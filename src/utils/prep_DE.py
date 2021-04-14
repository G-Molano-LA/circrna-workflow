#! /usr/bin/python3

# conda install numpy

import logging
from collection import namedtuple

LOGGER=logging.getLogger('prep_DE')
CIRC=namedtuple('CIRC','bsj fsj ratio rnaser_bsj rnaser_fsj')
INFO = namedtuple('INFO', 'strand circ_type gene_id gene_name gene_type')
"""
Python supports a type of container like dictionaries called “namedtuple()”
present in module, “collections“. Like dictionaries they contain keys that are
hashed to a particular value. But on contrary, it supports both access from key
value and iteration, the functionality that dictionaries lack.
"""
    LOGGER.info('Loading CIRIquant result: {}'.format(in_file))

    circ_data = {}
    circ_info = {}
    with open(in_file, 'r') as f:
        header = {}
        for line in f:
            if line.startswith('##'):
                key, value = line.rstrip().strip('#').split(':')
                header.update({key.strip(): value.strip()})
                continue

            content = line.rstrip().split('\t')
            tmp_parser = GTFParser(content)
            circ_data[tmp_parser.attr['circ_id']] = CIRC(
                float(tmp_parser.attr['bsj']),
                float(tmp_parser.attr['fsj']),
                float(tmp_parser.attr['junc_ratio']),
                float(tmp_parser.attr['rnaser_bsj']) if 'rnaser_bsj' in tmp_parser.attr else None,
                float(tmp_parser.attr['rnaser_fsj']) if 'rnaser_fsj' in tmp_parser.attr else None,
            )
            circ_info[tmp_parser.attr['circ_id']] = INFO(
                tmp_parser.strand,
                tmp_parser.attr['circ_type'],
                tmp_parser.attr['gene_id'] if 'gene_id' in tmp_parser.attr else 'NA',
                tmp_parser.attr['gene_name'] if 'gene_name' in tmp_parser.attr else 'NA',
                tmp_parser.attr['gene_type'] if 'gene_type' in tmp_parser.attr else 'NA',
            )
    return header, circ_data, circ_info


class GTFParser(object):
    """
    Class for parsing annotation gtf
    """

    def __init__(self, content):
        self.chrom = content[0]
        self.source = content[1]
        self.type = content[2]
        self.start, self.end = int(content[3]), int(content[4])
        self.strand = content[6]
        self.attr_string = content[8]


    @property
    def attr(self):
        """
        Parsing attribute column in gtf file
        """
        field = {}
        for attr_values in [re.split(r'[\s=]+', i.strip()) for i in self.attr_string.split(';')[:-1]]:
            key, value = attr_values[0], attr_values[1:]
            field[key] = ' '.join(value).strip('"')
        return field
