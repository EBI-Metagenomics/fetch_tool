#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import os
import subprocess
import requests

__author__ = "Hubert Denise, Maxim Scheremetjew"
__copyright__ = "Copyright (c) 2018 EMBL - European Bioinformatics Institute"
__version__ = "1.0"
__status__ = "Production"

"""
    Utility script which allows performing POST requests to ENA's rest api to fetch metadata for the various
    object/result types (e.g. sample, read_run or analysis)
"""


class InvalidFieldSuppliedException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MismatchedTokenException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class InvalidUsernamePasswordException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def get_sample_search_fields():
    return [
        "altitude",
        "broker_name",
        "center_name",
        "checklist",
        "collected_by",
        "collection_date",
        "country",
        "depth",
        "elevation",
        "environment_biome",
        "environment_feature",
        "environment_material",
        "environmental_package",
        "experimental_factor",
        "first_public",
        "host_body_site",
        "host_common_name",
        "host_genotype",
        "host_growth_conditions",
        "host_phenotype",
        "host_scientific_name",
        "host_sex",
        "host_status",
        "host_tax_id",
        "investigation_type",
        "last_updated",
        "location",
        "ph",
        "project_name",
        "region",
        "salinity",
        "sample_accession",
        "sample_alias",
        "sample_collection",
        "description",
        "secondary_sample_accession",
        "sequencing_method",
        "status_id",
        "target_gene",
        "temperature",
        "tax_id",
        "scientific_name",
    ]


def get_read_run_search_fields():
    return [
        "first_public",
        "instrument_model",
        "instrument_platform",
        "sequencing_method",
    ]


def get_analysis_search_fields():
    return ["first_public", "sequencing_method"]


def create_string_query(field_list):
    """
    generic function to transform a list of fields to return, separated by ", ", into a field query for ENA API
    :param field_list: list of accessions to be queried
    :return a string that can be directly integrated into a curl query.
    """
    if len(field_list) > 1:
        return_string = "%2C".join(field_list)
    else:
        return_string = field_list[0]
    return return_string


def create_sample_search_batch_query(secondary_sample_accs):
    """

    :param secondary_sample_accs: List of secondary sample accessions.
    :type secondary_sample_accs: List of str.
    """
    return (
        "secondary_sample_accession%3D"
        + "%20OR%20secondary_sample_accession%3D".join(secondary_sample_accs)
    )


def create_study_search_batch_query(secondary_study_accession):
    """

    :param secondary_study_accession: List of secondary study accessions.
    :type secondary_study_accession: List of str.
    """
    return "secondary_study_accession%3D" + "%20OR%20secondary_study_accession%3D".join(
        secondary_study_accession
    )


def perform_api_request(search_query):
    """

    :param search_query: Curl GET or POST request represented as a comma separated list.
    :type search_query: list
    :return: response body.
    @:rtype: bytes
    """
    with open(os.devnull, "w") as SNULL:
        try:
            raw_data_from_api = subprocess.check_output(search_query, stderr=SNULL)
        except:
            logging.error(
                "Could not get any sample metadata from ENA API. Program will exit now."
            )
            raise
    SNULL.close()

    return raw_data_from_api


def retrieve_metadata(result_type, query, search_fields, user_pass, api_url):
    """

    :param result_type: Check for valid result types ENA's api documentation. Valid values are for instance samples, read_run or analysis
    :param query: Search query, that could be multiple accessions or single accessions.
    :param api_url:
    :param user_pass:
    :param search_fields:
    :type result_type: str
    :type query: str
    :type search_fields: list
    :type api_url: str
    :type user_pass: str
    :return:
    """
    # transform search fields
    search_fields_str = create_string_query(search_fields)
    user_pass_chunks = user_pass.split(":")
    username = user_pass_chunks[0] if len(user_pass_chunks) == 2 else ""
    password = user_pass_chunks[1] if len(user_pass_chunks) == 2 else ""
    # create API curl query
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    headers = {
        "Accept": "text/plain",
        "Content-Type": "application/x-www-form-urlencoded",
    }
    data = {
        "dataPortal": "metagenome",
        "result": result_type,
        "query": query,
        "fields": search_fields_str,
    }

    r = requests.post(url, headers=headers, data=data, auth=(username, password))
    response_str = r.text
    if response_str.startswith("Invalid field(s) supplied"):
        raise InvalidFieldSuppliedException(response_str)
    elif response_str.startswith("MismatchedTokenException"):
        raise MismatchedTokenException(response_str)
    elif "Username/Password are invalid" in response_str:
        raise InvalidUsernamePasswordException(response_str)
    return response_str


def parse_response_str(response_str, accession_field):
    results = {}
    raw_lines = response_str.split("\n")
    # define the header line with the attributes
    header_line = raw_lines[0].split("\t")
    # go through the line of entries (len(raw_lines) -1 as the last line is empty)
    for i in range(1, len(raw_lines) - 1):
        attribute_value_dict = {}
        # separate the entries per attributes
        value_list = raw_lines[i].split("\t")
        # extract the attributes and their values and fill a dictionary for each sample
        for index, item in enumerate(header_line):
            attribute_value_dict[item.upper()] = value_list[index]
        accession = attribute_value_dict.get(accession_field)
        results[accession] = attribute_value_dict
    return results
