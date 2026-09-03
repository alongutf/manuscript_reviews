"""Fetch the sigma-38 (RpoS) sigmulon gene list from RegulonDB.

Writes metadata/regulondb_sigma38_regulon.txt (one gene name per line, with a
provenance header recording the RegulonDB release and retrieval date).

Run from the repo root or from scripts/:
    python scripts/fetch_regulondb_sigma38.py

Note on TLS: regulondb.ccg.unam.mx serves a valid GlobalSign OV certificate but
does not send the intermediate CA, so default verification fails. Rather than
disabling verification, this script fetches the intermediate from the URL named
in the leaf certificate's Authority Information Access extension and appends it
to the certifi bundle, so the chain is still verified up to a trusted root.
"""

import datetime
import json
import os
import re
import socket
import ssl
import subprocess
import tempfile
import urllib.request

API = "https://regulondb.ccg.unam.mx/graphql"
HOST = "regulondb.ccg.unam.mx"


def _ca_bundle():
    """Return a CA bundle path that can verify HOST, patching in the missing
    intermediate via AIA if the server does not supply it."""
    import certifi

    base = certifi.where()
    ctx = ssl.create_default_context(cafile=base)
    try:
        with ctx.wrap_socket(socket.create_connection((HOST, 443), timeout=30),
                             server_hostname=HOST):
            return base  # chain is complete, nothing to do
    except ssl.SSLCertVerificationError:
        pass

    # Grab the leaf, read its AIA "CA Issuers" URI, download that intermediate.
    unverified = ssl._create_unverified_context()
    with unverified.wrap_socket(socket.create_connection((HOST, 443), timeout=30),
                                server_hostname=HOST) as s:
        leaf_der = s.getpeercert(True)

    tmp = tempfile.mkdtemp(prefix="regulondb_tls_")
    leaf = os.path.join(tmp, "leaf.der")
    with open(leaf, "wb") as fh:
        fh.write(leaf_der)
    text = subprocess.run(
        ["openssl", "x509", "-in", leaf, "-inform", "DER", "-noout",
         "-ext", "authorityInfoAccess"],
        capture_output=True, text=True, check=True).stdout
    m = re.search(r"CA Issuers - URI:(\S+)", text)
    if not m:
        raise RuntimeError("no CA Issuers URI in the server certificate")
    inter_der = urllib.request.urlopen(m.group(1), timeout=60).read()
    inter = os.path.join(tmp, "inter.der")
    with open(inter, "wb") as fh:
        fh.write(inter_der)
    inter_pem = os.path.join(tmp, "inter.pem")
    subprocess.run(["openssl", "x509", "-in", inter, "-inform", "DER",
                    "-out", inter_pem], check=True)

    bundle = os.path.join(tmp, "bundle.pem")
    with open(bundle, "w") as out:
        for part in (base, inter_pem):
            with open(part) as fh:
                out.write(fh.read())

    # Fail loudly if the patched bundle still does not verify the real chain.
    ctx = ssl.create_default_context(cafile=bundle)
    with ctx.wrap_socket(socket.create_connection((HOST, 443), timeout=30),
                         server_hostname=HOST):
        pass
    return bundle


def gql(query, ctx):
    """POST a GraphQL query string to the RegulonDB API and return its "data" object.

    Raises if the response carries a top-level "errors" list (truncated to 500 chars
    so a malformed query doesn't dump a huge error payload to the console).
    """
    req = urllib.request.Request(
        API, data=json.dumps({"query": query}).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    out = json.loads(urllib.request.urlopen(req, context=ctx, timeout=300).read())
    if "errors" in out:
        raise RuntimeError(json.dumps(out["errors"])[:500])
    return out["data"]


def main():
    """Fetch the sigma-38 sigmulon and write the gene list with a provenance header.

    Resolves metadata/regulondb_sigma38_regulon.txt relative to the repo root (one
    level up from this script), builds a TLS context using the AIA-patched CA bundle
    from _ca_bundle(), then issues two GraphQL queries: one for the current RegulonDB
    release metadata (for the provenance header) and one searching sigmulons for
    "RpoS" to get sigma-38's gene list.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    dest = os.path.join(root, "metadata", "regulondb_sigma38_regulon.txt")

    ctx = ssl.create_default_context(cafile=_ca_bundle())

    releases = gql("{getDatabaseInfo{regulonDBVersion ecocycVersion releaseDate "
                   "genomeVersion}}", ctx)["getDatabaseInfo"]
    rel = releases[0]

    data = gql('{getSigmulonBy(search:"RpoS"){data{_id '
               'sigmaFactor{name abbreviatedName sigmulonGenes{_id name}} '
               'statistics{genes}}}}', ctx)["getSigmulonBy"]["data"]

    # the "RpoS" search can return other sigma factors/sigmulons that merely mention
    # RpoS, so pick out the record whose abbreviated name is exactly sigma38
    rec = next(r for r in data if r["sigmaFactor"]["abbreviatedName"] == "sigma38")
    sf = rec["sigmaFactor"]
    # de-duplicate and sort for a stable, diffable output file across reruns
    genes = sorted({g["name"] for g in sf["sigmulonGenes"]})
    if len(genes) != rec["statistics"]["genes"]:
        print("note: %d unique names vs statistics.genes=%d"
              % (len(genes), rec["statistics"]["genes"]))

    with open(dest, "w") as fh:
        fh.write("# sigma-38 (RpoS) sigmulon gene list\n")
        fh.write("# source: RegulonDB GraphQL API %s\n" % API)
        fh.write('# query: getSigmulonBy(search:"RpoS") -> sigmaFactor.sigmulonGenes\n')
        fh.write("# sigmulon id: %s (%s)\n" % (rec["_id"], sf["name"]))
        fh.write("# RegulonDB release: %s (released %s, EcoCyc %s, genome %s)\n"
                 % (rel["regulonDBVersion"], rel["releaseDate"],
                    rel["ecocycVersion"], rel["genomeVersion"]))
        fh.write("# retrieved: %s\n" % datetime.date.today().isoformat())
        fh.write("# n genes: %d\n" % len(genes))
        for g in genes:
            fh.write(g + "\n")

    print("RegulonDB %s (%s); wrote %d genes -> %s"
          % (rel["regulonDBVersion"], rel["releaseDate"], len(genes), dest))


if __name__ == "__main__":
    main()
