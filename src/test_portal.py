import requests
import os
import logging
import ftplib
import sys
from src.ENADAO import ENADAO


'''params = (
    ('dataPortal', 'metagenome'),
    ('dccDataOnly', 'false'),
    ('result', 'analysis'),
    ('format', 'tsv'),
    ('query', 'study_accession="PRJEB38205"'),
    ('fields', 'analysis_accession,study_accession,secondary_study_accession,sample_accession,'
               'secondary_sample_accession,analysis_title,analysis_type,center_name,first_public,last_updated,'
               'study_title,analysis_alias,study_alias,submitted_md5,submitted_ftp,generated_md5,generated_ftp,'
               'sample_alias,broker_name,sample_title'),
    ('download', 'true'),
)'''

'''def build_and_post_authenticated_search_request(self,
                                                result, query, fields,
                                                data_portal="metagenome",
                                                dcc_data_only=False,
                                                output_format="json",):
    """
        Builds and posts an authenticated search request.
    :return: server response
    """
    headers = {'Accept': '*/*', 'Content-Type': 'application/x-www-form-urlencoded'}
    data = {
        'dataPortal': data_portal,
        'dccDataOnly': dcc_data_only,
        'result': result,
        'query': query,
        'fields': fields,
        'format': output_format
    }
    r = requests.post(self.config['enaAPIUrl'] + 'search', headers=headers, data=data,
                      auth=(self.config['enaAPIUsername'], self.config['enaAPIPassword']))
    if r.status_code != 200:
        if r.status_code == 204:
            logging.warning("Could not retrieve any assemblies!")
        elif r.status_code == 401:
            logging.warning("Invalid Username or Password!")
        else:
            logging.warning("Received the following unknown response code from the "
                          "Portal API server:\n{}".format(r.status_code))
        logging.warning(self.PROGRAM_EXIT_MSG)
        sys.exit(1)
    else:
        return r'''


user = os.getenv('ENA_API_USER')
password = os.getenv('ENA_API_PASSWORD')
#NB. Original query string below. It seems impossible to parse and
#reproduce query strings 100% accurately so the one below is given
#in case the reproduced version is not "correct".

ENA_PORTAL_API_URL = 'https://www.ebi.ac.uk/ena/portal/api/search?dataPortal=metagenome&dccDataOnly=false&result='\
                     'analysis&query=study_accession=%22{0}%22&fields=analysis_accession,'\
                     'study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,'\
                     'analysis_title,analysis_type,center_name,first_public,last_updated,study_title,analysis_alias,'\
                     'study_alias,submitted_md5,submitted_ftp,generated_md5,generated_ftp,sample_alias,broker_name,'\
                     'sample_title&download=true'

ENA_OLD_ASSEMBLY_API_URL = 'https://www.ebi.ac.uk/ena/portal/api/search?result=wgs_set&query=study_accession=' \
                           '%22{0}%22&fields=study_accession,assembly_type,scientific_name,fasta_file&format=tsv'

project_acc = 'PRJEB22493'
#URL = ENA_PORTAL_FTP_URL.format(project_acc)

#response = requests.get(ENA_PORTAL_FTP_URL.format(project_acc), auth=(user, password))


def _retrieve_ena_url(url):
    attempt = 0
    response = None
    while True:
        try:
            # response = urllib.request.urlopen(url)
            response = requests.get(url, auth=(user, password))
            break
        except requests.exceptions.RequestException as e:  # check syntax
            logging.error(e)
            logging.warning("Error opening url " + url)
            attempt += 1
        if attempt >= 0:
            logging.critical("Failed to open url " + url + " after " + str(
                attempt) + " attempts")
            sys.exit(1)
    data = response.text.split('\n')
    headers = data[0].split('\t')
    data = list(filter(None, data))[1:]
    data = [{k: v for k, v in zip(headers, d.split('\t'))} for d in data]  # Convert rows to list of dictionaries
    return data


def _is_rawdata_filetype(filename):
    return any(x in filename for x in ['.fa', '.fna', '.fasta', '.fq', 'fastq'])


def _filter_secondary_files(joined_file_names, md5s):
    file_names = joined_file_names.split(';')
    md5s = md5s.split(';')
    filename_md5s = zip(file_names, md5s)
    filtered_filename_md5s = [(f, md5) for f, md5 in filename_md5s]
    print(filtered_filename_md5s)
    #[(f, md5) for f, md5 in filename_md5s if _is_rawdata_filetype(f)]
    filtered_file_names, filtered_md5s = zip(*filtered_filename_md5s)
    return filtered_file_names, filtered_md5s


def _rename_raw_files(file_names, run_id):
    print(file_names)
    file_names = [f.lower() for f in file_names]
    if any([".fastq" in fn for fn in file_names]):
        filetype = ".fastq.gz"
    elif any(x in ";".join(file_names) for x in ['.fasta', '.fna', '.fa']):
        filetype = ".fasta.gz"
    else:
        raise ValueError("Unknown sequence file format: " + ",".join(file_names))
    if len(file_names) == 1:
        return [run_id + filetype]
    else:
        return [run_id + '_' + str(i + 1) + filetype for i, _ in enumerate(file_names)]


def _get_raw_filenames(filepaths, md5s, run_id, is_submitted_file):
    filepaths, md5s = _filter_secondary_files(filepaths, md5s)
    if is_submitted_file:
        file_names = _rename_raw_files(filepaths, run_id)
    else:
        file_names = [os.path.basename(f) for f in filepaths]
    return filepaths, file_names, md5s


def map_datafields_ftp_2_db(assemblydata):
    is_submitted_file = assemblydata['submitted_ftp'] is not ''
    assemblydata['STUDY_ID'] = assemblydata.pop('secondary_study_accession')
    assemblydata['SAMPLE_ID'] = assemblydata.pop('secondary_sample_accession')
    assemblydata['DATA_FILE_PATH'] = assemblydata.pop('generated_ftp')
    assemblydata['ANALYSIS_ID'] = assemblydata.pop('analysis_accession')
    assemblydata['DATA_FILE_PATH'], assemblydata['file'], assemblydata['MD5'] = _get_raw_filenames(
        assemblydata['DATA_FILE_PATH'],
        assemblydata.pop('generated_md5'),
        assemblydata['ANALYSIS_ID'],
        is_submitted_file)
    return assemblydata


def _retrieve_project_info_from_api(project_accession):
    data = _retrieve_ena_url(ENA_PORTAL_API_URL.format(project_accession))
    logging.info("Retrieved {count} assemblies for study {project_accession} from "
                 "the ENA Portal API.".format(count=len(data), project_accession=project_accession))
    data = [d for d in data if d['analysis_type'] == 'SEQUENCE_ASSEMBLY']
    return list(map(map_datafields_ftp_2_db, data))


url = 'ftp.dcc-private.ebi.ac.uk/vol1/sequence/ERZ147/ERZ1478572/contigs.fa.gz'
dest = 'ERZ1478572.fasta.gz'


def download_lftp(dest, url):
    server = 'ftp.dcc-private.ebi.ac.uk'
    path_list = url.split('ebi.ac.uk/')[-1].split('/')[:-1]
    path = '/'.join(path_list)
    file_name = url.split('/')[-1]
    attempt = 0
    while attempt <= 3:
        try:
            with ftplib.FTP(server) as ftp:
                print("Downloading file from FTP server..." + url)
                print('Logging in...')
                ftp.login(user, password)
                ftp.cwd(path)
                print('Getting the file...')
                # store with the same name
                with open(dest, 'wb') as output_file:
                    ftp.retrbinary('RETR ' + file_name, output_file.write)
                print('File ' + dest + ' downloaded.')
                return True
        except ftplib.all_errors as e:
            print(e)
            attempt += 1
    else:
        print('Too many failed attempts. Trying wget now...')
        return False


'''if 'Login incorrect' in e:
    print('FTP login error')
elif 'Failed to change directory' in e:
    attempt += 1
    print('dir does not exist:' + path)
elif 'Failed to open file' in e:
    print('file does not exit:' + url)
    attempt += 1
if attempt >= self.config['url_max_attempts']:
    logging.critical("Failed to retrieve" + url + " after " + str(
        attempt) + " attempts")
    if self.interactive_mode:
        var = input(
            "Please type C to continue to fetch the next sequence file or anything else to exit: ")
        if not var.upper().startswith('C'):
            logging.info("Exiting now")
            sys.exit(0)
        else:
            break
    else:
        if self.force_mode:
            logging.warning(
                "Force mode is activated. Will skip the download of this run and move onto the next sequence!")
            break
        else:
            logging.warning(
                "Too many failed attempts. Program will exit now. " +
                "Try again to fetch the data in interactive mode (-i option)!")
            sys.exit(1)


def download_ftp(self, dest, url):
    if url[:4] == 'ftp.':
        url = 'ftp://' + url
    attempt = 0
    while True:
        try:
            logging.info("Downloading file from FTP server..." + url)
            download_command = ["wget", "-v" if self.args.verbose else "--user={}".format(ENA_API_USER),
                                "--password={}".format(ENA_API_PASSWORD), "-q", "-t", "5", "-O", dest, url]
            retcode = call(download_command)
            if retcode:
                logging.error("Error downloading the file from " + url)
            else:
                logging.info("Done.")
            break
        except IOError as err:
            logging.error("Error downloading the file from " + url)
            logging.error(err)
            attempt += 1
        if attempt >= self.config['url_max_attempts']:
            logging.critical("Failed to retrieve" + url + " after " + str(
                attempt) + " attempts")
            if self.interactive_mode:
                var = input(
                    "Please type C to continue to fetch the next sequence file or anything else to exit: ")
                if not var.upper().startswith('C'):
                    logging.info("Exiting now")
                    sys.exit(0)
                else:
                    break
            else:
                if self.force_mode:
                    logging.warning(
                        "Force mode is activated. Will skip the download of this run and move onto the next sequence!")
                    break
                else:
                    logging.warning(
                        "Too many failed attempts. Program will exit now. " +
                        "Try again to fetch the data in interactive mode (-i option)!")
                    sys.exit(1)

#ftplib.error_perm: 550 Failed to change directory.
#ftplib.error_perm: 530 Login incorrect.
#ftplib.error_perm: 550 Failed to open file.'''

#download_lftp(dest, url)


def _get_study_wgs_analyses(primary_project_accession):
    return ENADAO.retrieve_assembly_data(primary_project_accession)


def retrieve_insdc_assemblies(study_accession):
    # Step 1: Retrieve analysis results including assembly type from ENA Portal API
    json_data = _retrieve_ena_url(ENA_OLD_ASSEMBLY_API_URL.format(study_accession))
    print("Retrieved {count} assemblies of type wgs_set from ENA Portal API.".format(count=len(json_data)))
    if len(json_data):
        study_assembly_data = _get_study_wgs_analyses(study_accession)
        study_analyses = _combine_analyses(json_data, study_assembly_data)


@staticmethod
def _combine_analyses(metadata, wgs):
    combine_analyses = []
    wgs_dict = {}
    for assembly in wgs:
        analysis_id = assembly['ASSEMBLY_ID']
        wgs_dict[analysis_id] = {'CONTIG_CNT': assembly['CONTIG_CNT'],
                                 'FILENAME': assembly['WGS_ACC'],
                                 'GC_ID': assembly['GC_ID']}
    for analysis in metadata:
        analysis_id = analysis['accession']
        if analysis_id in wgs_dict:
            new_analysis = analysis.copy()
            new_analysis.update(wgs_dict[analysis_id])
            new_analysis['DATA_FILE_PATH'] = new_analysis.pop('fasta_file')
            combine_analyses.append(new_analysis)
    return combine_analyses


    # Convert JSON into a more handle format
    #study_analyses, status_ids = self.convert_api_response(json_data)
    #is_public = self.evaluate_statues(status_ids)
    #return self.retrieve_download_paths(study_accession, study_analyses, is_public)


    #return r.json()

#https://www.ebi.ac.uk/ena/portal/api/search?result=wgs_set&query=study_accession=%22PRJEB22493%22&fields=study_accession,assembly_type,scientific_name,fasta_file&format=tsv
# curl -X POST -H "Content-Type: application/x-www-form-urlencoded" -d 'result=wgs_set&query=study_accession=%22PRJEB22493%22&fields=study_accession,assembly_type,scientific_name,fasta_file&format=tsv' "https://www.ebi.ac.uk/ena/portal/api/search"


#print(_retrieve_ena_url(ENA_OLD_ASSEMBLY_API_URL))
retrieve_insdc_assemblies(project_acc)