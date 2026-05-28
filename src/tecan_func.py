import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import re
from scipy.ndimage import gaussian_filter1d

class tecan:
    """
    Class for handling Tecan data files.
    """

    def __init__(self, path_to_data, path_to_info):
        """
        Initializes the Tecan class with the paths to the data and info files.
        Args:
            path_to_data (str): Path to the data file.
            path_to_info (str): Path to the info file.
        """
        plate_info = pd.read_csv(path_to_info, index_col=0, header=0)
        self.plate_info = plate_info
        data = pd.read_excel(path_to_data)
        n_cycles = data.iloc[get_string_loc(data, 'Kinetic Cycles')[0][0], 4]
        channel_pos = get_string_loc(data, 'Cycle Nr.')
        channel_names = [data.iloc[pos[0] - 1, 0] for pos in channel_pos]
        self.channels = {name: data.iloc[pos[0]:pos[0]+n_cycles,:] for name, pos in zip(channel_names, channel_pos)}
        for name, df in self.channels.items():
            df.set_index(df.columns[0], inplace=True)
            df.columns = df.iloc[0]
            # remove nan
            df.dropna(axis=0, how='all', inplace=True)
            # replace non-numeric values with the previous value in the column
            self.channels[name] = replace_non_numeric(df[1:])

    def condition_locations(self, condition):
        """
        Extracts unique conditions and their corresponding well locations from a table of plate layout.

        Args:
            condition (str): name of condition to find the plate locations.

        Returns:
            locations (list): A list of strings representing the well location of the condition.
        """
        locations = []
        matches = self.plate_info.isin([condition])  # Boolean DataFrame
        locations = matches.stack()[matches.stack()]  # Keep only True values
        return [f"{row}{col}" for row, col in locations.index]

    def plot_condition(self,channel,condition):
        """
        Plot data from `channel` for one or multiple conditions.
        If multiple conditions are given, lines are colored by condition,
        and each location remains a separate line.

        Args:
            channel (str): Name of the channel to plot.
            condition (str | list[str]): Condition or list of conditions.

        Returns:
            fig: Plotly figure object.
        """
        # Normalize to list
        if isinstance(condition, str):
            conditions = [condition]
        else:
            conditions = list(condition)

        data = self.channels[channel]

        # Build a long-form DF: Time [s], value, location, condition, trace_id
        parts = []
        for cond in conditions:
            locs = self.condition_locations(cond)  # list of column names for this condition
            if not locs:
                continue  # skip empty condition gracefully
            sub = data[['Time [s]'] + locs].melt(
                id_vars='Time [s]',
                var_name='location',
                value_name='value'
            )
            sub['condition'] = cond
            # unique line per (condition, location)
            sub['trace_id'] = sub['condition'] + ' | ' + sub['location']
            parts.append(sub)

        if not parts:
            raise ValueError("No locations found for the provided condition(s).")

        df_long = pd.concat(parts, ignore_index=True)

        # Color by condition; separate line per location within each condition
        fig = px.line(
            df_long,
            x='Time [s]',
            y='value',
            color='condition',
            line_group='trace_id',
            hover_data={'location': True, 'condition': True, 'trace_id': False},
            title=", ".join(conditions) if len(conditions) > 1 else conditions[0]
        )
        fig.update_layout(yaxis_title=channel, legend_title_text='Condition')

        return fig


    def remove_background(self, fixed_background=0):
        """
        Removes the background optical density (OD) from OD channel in the Tecan data.
        This method assumes that the background data is in wells with condition=='MED'
        stores the fixed OD in a new channel: "OD_corrected"
        """
        OD_corrected = self.channels['OD'].copy()
        if fixed_background != 0:
            for col in OD_corrected.columns:
                if re.fullmatch(r"[A-Z]\d+", col):
                    OD_corrected[col] = OD_corrected[col] - fixed_background
                else:
                    continue
            self.channels['OD_corrected'] = OD_corrected
        else:
            background_loc = self.condition_locations('MED')
            background_avg = np.zeros(len(self.channels['OD']['Time [s]']))
            n = 0
            for well in background_loc:
                background_data = self.channels['OD'].loc[:, well]
                if np.mean(background_data[-100:])>0.1:
                    print(f"Warning: Background data for {well} is too high, skipping this well.")
                    continue
                else:
                    background_avg += background_data
                    n += 1
            for col in OD_corrected.columns:
                if re.fullmatch(r"[A-Z]\d+", col):
                    OD_corrected[col] = OD_corrected[col] - background_avg / n
                else:
                    continue

            self.channels['OD_corrected'] = OD_corrected


    def remove_spikes(self, channel, threshold=0.01):
        """
        Removes spikes from a specified channel in the Tecan data.
        A spike is defined as a change greater than the threshold between consecutive time points.
        The spike is replaced with the average of the two surrounding points.

        Args:
            channel (str): Name of the channel to process.
            threshold (float): Threshold for detecting spikes. Default is 0.01.
        """
        data = self.channels[channel].copy()
        for col in data.columns:
            if re.fullmatch(r"[A-Z]\d+", col):
                for i in range(2, len(data)-2):
                    if data.iloc[i, data.columns.get_loc(col)] - data.iloc[i-1, data.columns.get_loc(col)] > threshold and data.iloc[i+1, data.columns.get_loc(col)] - data.iloc[i, data.columns.get_loc(col)] < -threshold:
                        # replace spike with average of surrounding points
                        data.iloc[i, data.columns.get_loc(col)] = (data.iloc[i-1, data.columns.get_loc(col)] + data.iloc[i+1, data.columns.get_loc(col)]) / 2
        self.channels[channel] = data


    def align_times(self, channel, treatment_times):
        """
        Aligns the time points of a specified channel in the Tecan data based on treatment times.
        The time points are adjusted so that the treatment time is set to zero.

        Args:
            channel (str): Name of the channel to process.
            treatment_times (dict): Dictionary mapping conditions to their treatment times.
        """
        data = self.channels[channel].copy()
        values = self.plate_info.values.flatten()
        unique_conditions = pd.unique(values[~pd.isna(values)])
        max_ind = 0
        for times in treatment_times.keys():
            ind = np.argmin(np.abs(data['Time [s]'] - treatment_times[times]))+1
            if ind > max_ind:
                max_ind = ind
            for condition in unique_conditions:
                if times in condition:
                    for well in self.condition_locations(condition):
                        data[well] = data[well].shift(-ind)
                        data[well] = data[well]/data[well].iloc[0]
                else:
                    continue
        # number of rows to keep
        n = len(data)-max_ind
        self.channels[f"{channel}_aligned"] = data.iloc[:n, :]

    def get_derivative(self, channel, sigma=2):
        """
        Calculates the derivative of a specified channel in the Tecan data.
        The derivative is calculated as the difference between consecutive time points divided by the time difference.

        Args:
            channel (str): Name of the channel to process.
        """
        data = self.channels[channel].copy()
        derivative = data.copy()
        time = data['Time [s]'].to_numpy(dtype=float)
        for col in data.columns:
            if re.fullmatch(r"[A-Z]\d+", col):
                y = data[col].to_numpy(dtype=float)
                y_smoothed = gaussian_filter1d(y, sigma=sigma)
                dy_dx = np.gradient(y_smoothed, time)
                derivative[col] = dy_dx
            else:
                continue
        self.channels[f"{channel}_derivative"] = derivative

    def fluorescence_correction(self, channel):
        """
        Removes the background fluorescence from fluorescence channels in the Tecan data.
        stores the fixed fluorescence in a new channel: "{channel}_fixed"
        """
        values = self.plate_info.values.flatten()
        unique_conditions = pd.unique(values[~pd.isna(values)])
        fluorescence_corrected = self.channels[channel].copy()
        for condition in unique_conditions:
            if 'CTRL_T' in condition or 'TEST_T' in condition:
                strain, time, rep = condition.split('_')
                autofluorescence_data = self.channels[channel][self.condition_locations(f'{strain}_{rep}')]
                autofluorescence_avg = autofluorescence_data.mean(axis=1)
                #medium fluorescence
                medium_fluorescence = self.channels[channel][self.condition_locations(f'MED_{time}')]
                medium_fluorescence_avg = medium_fluorescence.mean(axis=1)
                #medium background
                medium_background = self.channels[channel][self.condition_locations(f'MED')]
                medium_background_avg = medium_background.mean(axis=1)
                # subtract the background and autofluorescence from the data
                for well in self.condition_locations(condition):
                    fluorescence_corrected[well] = fluorescence_corrected[well] - autofluorescence_avg - medium_fluorescence_avg + medium_background_avg
            else:
                continue
        self.channels[f'{channel}_corrected'] = fluorescence_corrected


    def normalize_fluorescence(self, channel):
        """
        Normalizes the fluorescence channels in the Tecan data.
        This method assumes that OD_corrected channel is present and fluorescence correction has been applied
        """
        # check if OD_corrected channel is present
        if 'OD_corrected' not in self.channels.keys():
            raise ValueError("OD_corrected channel is not present. Please run remove_background() first.")
        if f'{channel}_corrected' not in self.channels.keys():
            raise ValueError(f"{channel}_corrected channel is not present. Please run fluorescence_correction() first.")

        normalized_fluorescence = self.channels[f'{channel}_corrected'].copy()
        OD_corrected = self.channels['OD_corrected'].copy()
        # add a small value to avoid division by zero
        OD_corrected[OD_corrected <= 0] = 0.01
        normalized_fluorescence.iloc[:,2:] = normalized_fluorescence.iloc[:,2:]/OD_corrected.iloc[:,2:]
        self.channels[f'{channel}_normalized'] = normalized_fluorescence


    def get_value(self, channel, condition, time):
        """
        Retrieves the value from a specific channel and condition at a given time point.

        Args:
            channel (str): Name of the channel to retrieve data from.
            condition (str): Condition to retrieve data for.
            time (float): Time point to get the value for.

        Returns:
            pd.Series: Series containing the values at the specified time point.
        """
        if channel not in self.channels:
            raise ValueError(f"Channel '{channel}' not found in the data.")
        # find the row corresponding to the time
        measurement_times = self.channels[channel]['Time [s]']
        time_ind = np.argmin(np.abs(measurement_times - time))
        values_at_time = self.channels[channel].iloc[time_ind, :]
        return values_at_time[self.condition_locations(condition)].to_list()


    def plot_aligned_times(self, channel, condition1, time1, condition2, time2):
        """
        plots the aligned values of two conditions according to times of treatment
        :param channel: channel to plot
        :param condition1: name of first condition
        :param time1: time of first condition treatment
        :param condition2: name of second condition
        :param time2: time of second condition treatment
        :return: plotly figure object
        """
        # get the indices of the timepoints
        measurement_times = self.channels[channel]['Time [s]']
        ind1 = np.argmin(np.abs(measurement_times - time1))
        ind2 = np.argmin(np.abs(measurement_times - time2))
        data = self.channels[channel]
        n = len(data)-np.maximum(ind1, ind2)
        locations1 = self.condition_locations(condition1)
        data_condition1 = self.channels[channel].iloc[ind1:ind1+n,:]
        # normalize by the first value
        data_condition1.iloc[:,1:] = data_condition1.iloc[:,1:] / data_condition1.iloc[1,1:]
        locations2 = self.condition_locations(condition2)
        data_condition2 = self.channels[channel].iloc[ind2:ind2+n,:]
        data_condition2.iloc[:,1:] = data_condition2.iloc[:,1:] / data_condition2.iloc[1,1:]
        fig = px.line(data_condition1, x='Time [s]', y=locations1, title=f"Aligned {condition1} vs {condition2}")
        fig.update_traces(line=dict(color="blue"))
        for col in locations2:
            fig.add_trace(
                go.Scatter(
                    x=data_condition1['Time [s]'],
                    y=data_condition2[col],  # <- 1-D Series
                    mode='lines',
                    name=f'{condition2} – {col}',
                    legendgroup=condition2,  # groups traces in legend
                    line=dict(dash='dot', color='red')  # visually separate conditions
                )
            )
        fig.update_layout(yaxis_title=channel)
        return fig

# helper functions
def get_string_loc(dataframe, string):
    """
    Finds the location of a string in a DataFrame.
    Args:
        dataframe (pd.DataFrame): DataFrame to search in.
        string (str): String to find.
    Returns:
        list: List of indices where the string is found.
    """
    x, y = np.where(dataframe == string)
    return list(zip(x, y))


def is_number(x):
    """
    Checks if the input can be converted to a number.
    :param x: Input value to check.
    :return: boolean
    """
    try:
        float(x)
        return True
    except (ValueError, TypeError):
        return False



def replace_non_numeric(df):
    """
    Replaces non-numeric values in a DataFrame with the previous value in the column.
    :param df: pandas dataframe to process.
    :return: df
    """
    for col in range(len(df.columns)):
        for i in range(len(df)):
            if not is_number(df.iloc[i, col]):
                if i == 0:  # first row
                    df.iloc[i, col] = 0
                else:
                    df.iloc[i, col] = df.iloc[i - 1, col]
    return df
