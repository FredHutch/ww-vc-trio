#https://developer.basespace.illumina.com/docs/content/documentation/native-apps/spacedock-conventions#UploadResultstoBaseSpacewithProperties

import glob, json, os, subprocess, sys

#most files are in /opt/files, so switch to it and use as a working directory
os.chdir('/opt/files')

print('Chromoseq Basespace App v1.0.0')
print('Source code available at https://github.com/genome/docker-basespace_chromoseq')

#download and setup VEP cache
print('Downloading VEP cache, this may take a while...')
subprocess.check_call(['curl', '-OsS', 'ftp://ftp.ensembl.org/pub/release-90/variation/VEP/homo_sapiens_vep_90_GRCh38.tar.gz'])
print('Download complete. Unpacking...')
subprocess.check_call(['tar', 'xzf', 'homo_sapiens_vep_90_GRCh38.tar.gz'])
print('Unpacking complete')


#get basespace related properties/metadata
with open("/data/input/AppSession.json") as a_s_j:
    appsession = json.load(a_s_j)

#finding the ID of the project from which this analysis was launched; this is needed in order
#to create the directory structure specified by basespace for automatic result uploading
for e in appsession['Properties']['Items']:
    if e['Name'] == 'Input.project-id':
        project_id = str(e['Content']['Id']) #note that this will return a unicode object, not a str; this is still python 2.7.x, so these are 2 different things
    elif e['Name'] == 'Input.ref-fa-id':
        ref_fa_href = str(e['Content']['Href'])
    elif e['Name'] == 'Input.app-result-id':
        dragen_href = str(e['Content']['Href'])
    elif e['Name'] == 'Input.gender-select-id':
        sample_sex = str(e['Content'])
appsession_href = appsession['Href'] #basespace internal reference to the current appsession

cram_search = glob.glob('/data/input/appresults/*/*.cram')
if len(cram_search) != 1:
    print('Error- expected 1 cram file but found {0}'.format(len(cram_search)))
    sys.exit(1)
crai_search = glob.glob('/data/input/appresults/*/*.crai')
if len(crai_search) != 1:
    print('Error- expected 1 crai file but found {0}'.format(len(cram_search)))
    sys.exit(1)
ref_search = glob.glob('/data/input/appresults/*/*.fa')
if len(ref_search) != 1:
    print('Error- expected 1 reference fasta file but found {0}'.format(len(ref_search)))
    sys.exit(1)

map_summary_search = glob.glob('/data/input/appresults/*/*.mapping_metrics.csv') #required
if len(map_summary_search) != 1:
    print('Error- expected 1 mapping summary file but found {0}'.format(len(map_summary_search)))
    sys.exit(1)

cov_summary_search = glob.glob('/data/input/appresults/*/*.wgs_coverage_metrics.csv') #optional
cov_found = len(cov_summary_search) == 1

cram_file = cram_search[0]
crai_file = crai_search[0]
ref_file_temp = ref_search[0]
name = cram_file.split("/")[-1].split(".")[0]

mapping_summary_file = map_summary_search[0]
if cov_found:
    coverage_summary_file = cov_summary_search[0]

ref_file = '/opt/files/all_sequences.fa'
subprocess.check_call(['ln', '-s', ref_file_temp, ref_file])

#basespace-specified directory structure: /data/output/appresults/<project-id>/[directory_with_appresult_name]/[your_files]
output_dir = "/data/output/appresults/{0}/{1}".format(project_id, name) 
os.makedirs(output_dir) #recursively create the output directory

#create inputs.json file
wf_inputs_dict = \
{
    "ChromoSeq.Cram": cram_file,
    "ChromoSeq.CramIndex": crai_file,
    "ChromoSeq.Name": name,
    "ChromoSeq.OutputDir": output_dir,
    "ChromoSeq.Reference": ref_file,
    "ChromoSeq.Translocations": "/opt/files/chromoseq_translocations.bedpe",
    "ChromoSeq.GenesBed": "/opt/files/chromoseq_genes.bed",
    "ChromoSeq.Cytobands": "/opt/files/hg38.cytoBandIdeo.bed.gz",
    "ChromoSeq.SVDB": "/opt/files/chromoseq_sv_filter.bedpe.gz",
    "ChromoSeq.MantaConfig": "/opt/files/configManta.hg38.py.ini",
    "ChromoSeq.ReferenceIndex": "/opt/files/all_sequences.fa.fai",
    "ChromoSeq.ReferenceBED": "/opt/files/all_sequences.fa.bed.gz",
    "ChromoSeq.VEP": "/opt/files",
    "ChromoSeq.gcWig": "/usr/local/lib/R/site-library/ichorCNA/extdata/gc_hg38_500kb.wig",
    "ChromoSeq.mapWig": "/usr/local/lib/R/site-library/ichorCNA/extdata/map_hg38_500kb.wig",
    "ChromoSeq.ponRds": "/opt/files/nextera_hg38_500kb_median_normAutosome_median.rds_median.n9.gr.rds",
    "ChromoSeq.centromeres": "/usr/local/lib/R/site-library/ichorCNA/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
    "ChromoSeq.genomeStyle": "UCSC",
    "ChromoSeq.genome": "hg38",
    "ChromoSeq.tmp": "/tmp",
    "ChromoSeq.minVarFreq": 0.02,
    "ChromoSeq.JobGroup": "dummy",
    "ChromoSeq.chromoseq_docker": "dummy",
    "ChromoSeq.CustomAnnotationVcf": "/opt/files/chromoseq_custom_anntations.040920.vcf.gz",
    "ChromoSeq.CustomAnnotationIndex": "/opt/files/chromoseq_custom_anntations.040920.vcf.gz.tbi",
    "ChromoSeq.HotspotVCF": "/opt/files/chromoseq_hotspot.vcf.gz",
    "ChromoSeq.CustomAnnotationParameters": "MYELOSEQ,vcf,exact,0,TCGA_AC,MDS_AC,MYELOSEQBLACKLIST",
    "ChromoSeq.Gender": sample_sex,
    "ChromoSeq.MappingSummary": mapping_summary_file
}

if cov_found:
    wf_inputs_dict['ChromoSeq.CoverageSummary'] = coverage_summary_file

with open("/opt/files/inputs.json", "w+") as f:
    json.dump(wf_inputs_dict, f)

#create metadata file required by basespace for upload; each workflow generated ouput file is tagged
#with the downloaded basespace file(s) used to generate it

metadata_outfile = output_dir + "/_metadata.json"

#note that the value of key "Name" must match the dirname at the trailing end of $output_dir
metadata_json_template = \
{
    "Name": name,
    "Description": "Outputs from Chromoseq workflow run on {}".format(cram_file),
    "HrefAppSession": appsession_href,
    "Properties": [
        {
            "Type": "appresult[]",
            "Name": "Input.AppResults",
            "Items": [dragen_href]
        },
        {
            "Type": "file[]",
            "Name": "Input.Files",
            "Items": [ref_fa_href]
        }
    ]
}

with open(metadata_outfile, "w+") as f:
    json.dump(metadata_json_template, f)

print('\nLaunching cromwell')
cromwell_cmd = ["/usr/bin/java", "-Dconfig.file=/opt/files/basespace_cromwell.config", "-jar", "/opt/cromwell-36.jar", "run", "-t", "wdl", "-i", "/opt/files/inputs.json", "/opt/files/Chromoseq.v17.wdl"]
subprocess.check_call(cromwell_cmd)

print('\nCromwell complete')
