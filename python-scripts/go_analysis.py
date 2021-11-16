#!/usr/bin/python
#-*-coding: utf-8-*-


from go_analysis_parser import go_parser
import go_analysis_functions as go_func
import platform
import subprocess as subp
import time

args = go_parser()

if platform.system() == 'Windows':
    subp.call('cls', shell=True)
elif platform.system() == 'Linux':
    subp.call('clear', shell=True)

start_time = time.time()

if args.task == 'download_go':
    print("\t\t\t=== Download task ===")
    print('\n' + time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + "Start")
    go_func.download_go(args.outFile)
if args.task == 'parse_go':
    print("\t\t\t=== Parsing  task ===")
    print('\n' + time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + "Start")
    go_func.parse_go(args.file, args.outdir)
if args.task == 'get_annotations':
    print("\t\t\t=== Annotations task ===")
    print('\n' + time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + "Start")
    go_func.treat_one_request("ensemblgenomes",args.org_name, args.outdir)

#    go_func.get_annotations(args.org_name, args.outdir)
if args.task == 'expand':
    print("\t\t\t=== Expansion task ===")
    print('\n' + time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + "Start")
    go_func.expand(args.org_name, args.description_file, args.relation_file, args.annotation_file)
if args.task == 'enrichment':
    print("\t\t\t=== Enrichment task===")
    print('\n' + time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + "Start")
    go_func.enrichment()

print('\n' + time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + "Job done")
print("\n--- %.*f seconds ---" % (3, (time.time() - start_time)))
