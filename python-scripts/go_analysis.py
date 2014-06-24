#!/usr/bin/python
#-*-coding: utf-8-*-

from go_analysis_parser import go_parser
import go_analysis_functions as go_func

args = go_parser()

if args.task == 'download_go':
    go_func.download_go(args)
if args.task == 'parse_go':
    go_func.parse_go()
if args.task == 'get_annotations':
    go_func.get_annotations(args.org_name)
if args.task == 'expand':
    go_func.expand(args.org_name)
if args.task == 'enrichment':
    go_func.enrichment()
if args.task == 'ancestor':
    go_func.ancestor(args)