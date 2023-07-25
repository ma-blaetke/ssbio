from os import path as op

import requests
import json
from lxml import etree

import ssbio.utils
import logging
log = logging.getLogger(__name__)

def blast_pdb(seq, outfile='', outdir='', evalue=0.0001, seq_ident_cutoff=0.0, link=False, force_rerun=False):
    """Returns a list of BLAST hits of a sequence to available structures in the PDB.

    Args:
        seq (str): Your sequence, in string format
        outfile (str): Name of output file
        outdir (str, optional): Path to output directory. Default is the current directory.
        evalue (float, optional): Cutoff for the E-value - filters for significant hits. 0.001 is liberal, 0.0001 is stringent (default).
        seq_ident_cutoff (float, optional): Cutoff results based on percent coverage (in decimal form)
        link (bool, optional): Set to True if a link to the HTML results should be displayed
        force_rerun (bool, optional): If existing BLAST results should not be used, set to True. Default is False

    Returns:
        list: Rank ordered list of BLAST hits in dictionaries.

    """

    if len(seq) < 12:
        raise ValueError('Sequence must be at least 12 residues long.')
    if link:
        page = 'PDB results page: http://www.rcsb.org/pdb/rest/getBlastPDB1?sequence={}&eCutOff={}&maskLowComplexity=yes&matrix=BLOSUM62&outputFormat=HTML'.format(seq, evalue)
        print(page)

    #parser = etree.XMLParser(ns_clean=True)

    outfile = op.join(outdir, outfile)
    if ssbio.utils.force_rerun(force_rerun, outfile):
        # Load the BLAST XML results if force_rerun=True
        page = f"https://search.rcsb.org/rcsbsearch/v2/query?json=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22sequence%22%2C%22parameters%22%3A%7B%22evalue_cutoff%22%3A{evalue}%2C%22identity_cutoff%22%3A{seq_ident_cutoff}%2C%22sequence_type%22%3A%22protein%22%2C%22value%22%3A%22{seq}%22%7D%7D%2C%22request_options%22%3A%7B%22scoring_strategy%22%3A%22sequence%22%7D%2C%22return_type%22%3A%22polymer_instance%22%7D"
    
        req = requests.get(page)
        if req.status_code == 200:
            response = req.json()['result_set']

            # Save the XML file
            if outfile:
                with open(outfile, 'w') as f:
                    json.dump(response, f)

            # Parse the XML string
            #tree = etree.ElementTree(etree.fromstring(response, parser))
            #log.debug('Loaded BLAST results from REST server')
        else:
            log.error('BLAST request timed out')
            return []
    else:
        #tree = etree.parse(outfile, parser)
        with open(outfile, 'r') as f:
            response = json.load(f)

        log.debug('{}: Loaded existing BLAST results'.format(outfile))

    hit_list = []

    for hit in response:

        info = {}

        info['hit_pdb'] = hit['identifier']
        info['hit_pdb_chains'] = hit['identifier']
        info['hit_score'] = hit['score']

        hit_list.append(info)

    log.debug("{}: Number of BLAST hits".format(len(hit_list)))
    return hit_list