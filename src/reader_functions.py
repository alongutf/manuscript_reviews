import pandas as pd
import numpy as np
import datetime
import colorsys


class DistinctColorGenerator:
    """
    Generate distinct colors
    Args:
        None
    Returns:
        r,g,b (int): RGB values of the generated color
    """
    def __init__(self):
        self.index = 0  # Keep track of the current color index

    def get_new_color(self):
        # Golden ratio conjugate for even color spacing
        golden_ratio_conjugate = 0.618033988749895

        # Compute a new hue value
        hue = (self.index * golden_ratio_conjugate) % 1.0
        self.index += 1  # Increment index for the next call

        # Convert HSL to RGB (S=0.8, V=0.8 for vibrant colors)
        r, g, b = colorsys.hsv_to_rgb(hue, 0.8, 0.8)
        return r, g, b


def extract_condition_locations(file_path):
    """
    Extracts unique conditions and their corresponding well locations from a table of plate layout.

    Args:
        file_path (str): Path to the CSV file containing the table.

    Returns:
        dict: A dictionary where keys are unique conditions and values are lists of well locations.
    """
    # Load the data
    table = pd.read_csv(file_path, header=None)

    # Extract column labels (from the first row) and row labels (from the first column)
    table_cleaned = table.iloc[1:, 1:]  # Exclude the first row and column
    table_cleaned.columns = table.iloc[0, 1:]  # Set column headers from the first row
    table_cleaned.index = table.iloc[1:, 0]  # Set row labels from the first column

    # Dictionary to store conditions and their locations
    condition_locations = {}

    # Iterate over the table to collect well locations for each condition
    for row_label, row in table_cleaned.iterrows():
        for col_label, condition in row.items():
            if pd.notna(condition):  # Only consider non-NaN conditions
                location = f"{row_label}{col_label}"
                if condition not in condition_locations:
                    condition_locations[condition] = []
                condition_locations[condition].append(location)

    return condition_locations


def get_read_data(df, key, read=1):
    """
    Extracts data for a given read from the data file.
    Args:
        df (pandas.DataFrame): Data frame containing the data.
        key (str): Key that defines the location of the sub-table that contains the data ('Time': OD measurements,
         'Wavelength': spectrum).
        read (int): index of the read to extract data from.
    Returns:
        selected_df (numpy.array): Dataframe that contains only the required data.
    """
    # Find the locations (row, column) where the key appears
    locations = [(row, col) for row in range(df.shape[0]) for col in range(df.shape[1]) if key == str(df.iat[row, col])]
    read_location = locations[read-1]
    # find the end of the sub-table
    end_index = read_location[0]
    while end_index < df.shape[0] and not pd.isna(df.iat[end_index, read_location[1]]):
        end_index += 1
    # extract the sub-table
    selected_df = df.iloc[read_location[0]:end_index, read_location[1]:read_location[1] + 96]
    selected_df.set_index(selected_df.columns[0], inplace=True)
    # set the first row as the column names
    selected_df.columns = selected_df.iloc[0]
    selected_df = selected_df[1:]
    time = [convert_to_minutes(t) for t in selected_df.index.values]
    selected_df.index = time
    # remove rows with index = 0
    selected_df = selected_df[selected_df.index != 0]
    return selected_df


def time_to_timedelta(t):
    time_delta = datetime.timedelta(hours=t.hour, minutes=t.minute, seconds=t.second, microseconds=t.microsecond)
    #convert to minutes
    return time_delta.total_seconds() / 60

def concat_read_data(df1, df2, time1, time2):
    """
    Concatenates two dataframes with different time points.
    Args:
        df1 (pandas.DataFrame): First data frame.
        df2 (pandas.DataFrame): Second data frame.
        time_diff (float): Time difference between the two data frames.
    Returns:
        concatenated_df (pandas.DataFrame): Concatenated data frame.
    """
    # shift the time index of the second dataframe
    time_diff = time_to_timedelta(time2) - time_to_timedelta(time1)
    df2.index = df2.index + time_diff
    # concatenate the two dataframes
    concatenated_df = pd.concat([df1, df2])
    return concatenated_df


def get_tecan_data(df, channel):
    """
    Extracts data for a given read from the data file.
    Args:
        df (pandas.DataFrame): Data frame containing the data.
        channel (str): Key that defines the location of the sub-table that contains the data ('OD': OD measurements,
         'SytoxBlue': spectrum).
    Returns:
        selected_df (numpy.array): Dataframe that contains only the required data.
    """
    # Find the locations (row, column) where the key appears
    locations = [(row, col) for row in range(df.shape[0]) for col in range(df.shape[1]) if channel == str(df.iat[row, col])]
    read_location = locations[0]
    # find the end of the sub-table
    end_index = read_location[0]+1
    while end_index < df.shape[0] and not pd.isna(df.iat[end_index, read_location[1]+1]):
        end_index += 1
    # extract the sub-table
    selected_df = df.iloc[read_location[0]+1:end_index, read_location[1]:read_location[1] + 63]
    selected_df.set_index(selected_df.columns[0], inplace=True)
    # set the first row as the column names
    selected_df.columns = selected_df.iloc[0]
    selected_df = selected_df[1:]
    return selected_df


def get_spectrum_data(path_to_data, path_to_info, condition, read=1):
    """
    Extracts spectral data for a given condition from data file.
    Args:
        path_to_data (str): Path to the data file.
        path_to_info (str): Path to the info file.
        condition (str): Condition for which to extract the data.
        read (int): index of the read to extract data from.
    Returns:
        wave_length (numpy.array): Array of wave lengths.
        intensity (numpy.array): Array of intensity values.
    """
    # Load the data
    data = pd.read_excel(path_to_data)
    # extract the read data
    data = get_read_data(data, 'Wavelength', read)

    # get locations of the condition
    locations = extract_condition_locations(path_to_info)
    locations = locations[condition]
    # create numpy arrays to store wave lengths and intensity values
    wave_length = data.index.values
    intensity = np.zeros((len(wave_length), len(locations)))
    # Extract the wave lengths and intensity values
    for i, location in enumerate(locations):
        intensity[:, i] = data[location].to_numpy()

    return wave_length, intensity


def convert_to_minutes(time):
    # convert time object to minutes
    time = str(time)
    day = 0
    if 'day' in time:
        time = time.split(', ')[-1]
        hour, minute, second = time.split(':')
        day = 24
    else:
        hour, minute, second = time.split(':')

    return (float(hour)+day) * 60 + float(minute) + float(second) / 60


def get_od_data(path_to_data, path_to_info, condition, read=1):
    """
    Extracts OD data for a given condition from data file.
    Args:
        path_to_data (str): Path to the data file.
        path_to_info (str): Path to the info file.
        condition (str): Condition for which to extract the data.
    Returns:
        times (numpy.array): Array of time points.
        od (numpy.array): Array of OD values.
    """
    # Load the data
    od_data = pd.read_excel(path_to_data)
    # extract the read data
    od_data = get_read_data(od_data, 'Time', read+1)
    # get locations of the condition
    locations = extract_condition_locations(path_to_info)
    locations = locations[condition]
    # find number of time points
    end_index = 0
    times = [convert_to_minutes(time) for time in od_data.index]
    for time in times:
        if time == 0:
            break
        end_index += 1

    od = np.zeros((end_index, len(locations)))
    # Extract the time points and OD values
    for i, location in enumerate(locations):
        od[:, i] = od_data[location].to_numpy()[:end_index]

    return times[:end_index], od

def shift_medium_fluoresence(medium_data, exp_data):
    """
    Shifts the medium fluorescence to match the time point of the EXP_T0 data.
    :param medium_data: data of the medium fluorescence
    :param exp_data: data of the exponential cell fluorescence
    :return: medium_data: shifted medium fluorescence data
    """
    ind_t0 = np.where(exp_data.values >1000)[0][0]
    int_t2 = np.where(medium_data.values > 1000)[0][0]
    ind_diff = int_t2 - ind_t0
    return np.roll(medium_data, -ind_diff)
