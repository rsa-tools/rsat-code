#Begin by importing the Requests module
import requests, sys
 
#State the base URL
server = "http://rsat-tagc.univ-mrs.fr/rest.wsgi"
#The endpoint indicates which RSAT resource you are interested in
ext = "/convert-seq/" ##Converts sequences between different formats. Optionally, also returns the reverse-complement of the input sequences, or perform some cleaning operations (skip short sequences, suppress Ns, …).

#Write the parameters specifying details of how you want to interact with the resource. For default option write None
data =  {
        "i_string" : " ", ##Input string specifying the query. The value can be the query content, the URL of a file available on some Web server, the internal path of the result file returned by another tool of this RSAT server (piping for workflows).
        "i_string_type" : "url", ##Type of information provided by the input string. Supported values: url: URL (Web address) to the input file; piping: result file from other tool; text: input content
        "mask" : None, ##String. Masked characters are replaced by by N characters, or by a dot (option -dotmask). Supported. upper, lower, non-dna
        "noempty" : None, ##Boolean. Remove empty sequences from the set (same as -skip_short 1).
        "mask_short" : None, ##Integer. Mask (replace by N characters) sequences shorter than the specified length.
        "skip_short" : None, ##Integer. Skip sequences shorter than the specified length. Same functionality as -mask_short, except that short sequences are not returned at all in the output.
        "skip_long" : None, ##Integer.  Skip sequences longer than the specified length.
        "last" : None, ##Integer. Stop after the Nth sequence.
        "top" : None, ##Integer. Same as -last N
        "first" : None, ##Integer. Start at the Nth sequence (skip the N-1 first sequences).
        "skip" : None, ##Integer. Skip the N first sequences (start at sequence N+1).
        "from" : None, ##String. Input format. Supported - embl, fasta, filelist, ft, gcg, genbank, ig, maf,multi, ncbi, raw,tab,wc, wconsensus
        "to" : None, ##String. Output format. Supported - fasta, fastq, filelist, ft, ig, multi, raw, tab, wc, wconsensus
        "id_col" : None, ##Integer. Column containing sequence identifiers in tab format.
        "seq_col" : None, ##Integer. Column containing sequence sequences in tab format.
        "comment_col" : None, ##Integer. Column containing sequence comments (sequence description) in tab format.
        "lw" : None, ##Integer. Line width. A carriage return is inserted every lw characters within the output sequence.
        "addrc" : None ##Boolean. Adds the reverse complement of each input sequence to the output file.
        "lc" : None ##Boolean. Lowercase. the sequence is printed in lowercase.
        "uc" : None ##Boolean. Uppercase. the sequence is printed in uppercase.
        "dna" : None ##Boolean. Convert any non-acgt character into “n” characters.
        "dotmask" : None ##Boolean. Convert masked characters into dots.
        "id" : None ##String. Identifier. sequence identifier (useful for converting a raw sequence from the STDIN)
        "prefix" : None ##String. Sequence prefix (useful for converting from a multi sequence)
        "nocheckid" : None ##Boolean. Prevent to check sequence IDs for conversion to file list.
} 
r = requests.get(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"}) ##Default value : text/plain
#r = requests.post(server+ext, data, headers={ "Content-Type" : "text/plain", "Accept" : "application/json"})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()

 
print(r.text)
#print (repr(r.json))