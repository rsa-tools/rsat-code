#!/usr/bin/perl
#### this cgi script fills the HTML form for the publications page
if ($0 =~ /([^(\/)]+)$/) {
    push (@INC, "$`lib/");
}

use CGI;
use CGI::Carp qw/fatalsToBrowser/;

require "RSA.lib";
require "RSA2.cgi.lib";
$ENV{RSA_OUTPUT_CONTEXT} = "cgi";

$query = new CGI;

$default{year} = '2017';
$default{type} = 'all';

foreach $key (keys %default) {
    if ($query->param($key)) {
        $default{$key} = $query->param($key);
    }
}

### print the form ###
&RSA_header("Publications", 'info');

print '<center>
<h2>
<a href="RSAT_home.cgi">RSAT</A> - Publications</H2>
</center>
';

print '<br/><span style="color: #cc6600;"><i class="fa fa-pencil fa-lg"></i> <b>Citing RSAT complete suite of tools:</b></span>
<div align="left" style="border:1px solid #F6E6CA; padding:10px">
<ul>
<li>Medina-Rivera A, Defrance M, Sand O, Herrmann C, Castro-Mondragon J, Delerce J, Jaeger S, Blanchet C, Vincens P, Caron C, Staines DM, Contreras-Moreira B, Artufel M, Charbonnier-Khamvongsa L, Hernandez C, Thieffry D, Thomas-Chollier M, van Helden J (2015)
<b>RSAT 2015: Regulatory Sequence Analysis Tools </b>. Nucleic Acids Res. 2015 (Web Server issue) in press.
<a href="http://nar.oxfordjournals.org/content/early/2015/04/21/nar.gkv362.full">[Full text]</a>
</li>
<li>Thomas-Chollier M, Defrance M, Medina-Rivera A, Sand O, Herrmann C, Thieffry D, van Helden J. (2011)
<b>RSAT 2011: regulatory sequence analysis tools</b>. Nucleic Acids Res. 2011 Jul;39(Web Server issue):W86-91.
<a href=http://www.ncbi.nlm.nih.gov/pubmed/21715389>[Pubmed 21715389]</a>
<a href="http://nar.oxfordjournals.org/content/39/suppl_2/W86.long">[Full text]</a>
</li>
<li>Thomas-Chollier, M., Sand, O., Turatsinze, J. V., Janky, R.,
Defrance, M., Vervisch, E., Brohee, S. & van Helden, J. (2008). <b>RSAT:
regulatory sequence analysis tools</b>. Nucleic Acids
Res.
<a href=http://www.ncbi.nlm.nih.gov/pubmed/18495751>[Pubmed 18495751]</a>
<a href="http://nar.oxfordjournals.org/cgi/content/full/36/suppl_2/W119">[Full text]</a>
</li>
<li>van Helden, J. (2003). <b>Regulatory sequence analysis
tools</b>. Nucleic Acids Res. 2003 Jul 1;31(13):3593-6. [<a target=_blank
href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=PubMed&list_uids=12824373&dopt=Abstract">Pubmed
12824373</a>] [<a target=_blank
href=http://nar.oupjournals.org/cgi/content/full/31/13/3593>Full
text</a>] [<a target=_blank
href=http://nar.oupjournals.org/cgi/reprint/31/13/3593.pdf>pdf</a>]
</li>
</ul>
</div>
';

my $pub_file = "publications.csv";
open(my $data, '<', $pub_file) or die "Could not open '$file' $!\n";
my %years;
while( my $line = <$data> ){
    chomp $line;
    
    my @fields = split ";", $line;
    if($years{$fields[0]}){
        push @{$years{$fields[0]}}, [$fields[1] . '. ' . $fields[2], $fields[3]];
    }else{
        $years{$fields[0]} = ();
        push @{$years{$fields[0]}}, [$fields[1] . '. ' . $fields[2], $fields[3]];
    }
}

my @years = sort {$b <=> $a} keys %years;
unshift @years,('All');

print "<br/><button onclick='show(\"all\")'>View all articles</button>&nbsp;&nbsp; Year: ";
print $query->popup_menu(-name=>'year', -id=>'year', -Value=>[@years], -default=>$default{year}, -onchange=>'show()');
print "&nbsp;&nbsp;Type: ";
print $query->popup_menu(-name=>'type', -id=>'type', -values=>['All','RSAT suite', 'Individual tools','Protocols','Others'], -default=>'All', onchange=>'show()');

print "<div id='result' style='border:solid 1px #F6E6AC; padding:10px; margin-top:10px'>";
print "<h3>2017</h3>";
print "<ol>";
foreach $_ (@{$years{'2017'}}){
        print "<p><li>";
        print ${$_}[0];
        print "</li></p>";
}
print "</ol>";


print "<script type='text/javascript'>
function show(sel) {
    year = document.getElementById('year').value;
    type = document.getElementById('type').value;
    if(type == 'All'){ type = 'all';}
    else if (type == 'RSAT suite'){ type = 'suite';}
    else if (type == 'Individual tools'){ type = 'tools';}
    else if (type == 'Protocols'){ type = 'protocols';}
    else if (type == 'Others'){ type = 'others';}

if(sel == 'all'){
    year = 'All';
    type='all';
    document.getElementById('year').value = 'All';
    document.getElementById('type').value = 'All';
}

    \$.ajax({
        type: 'GET',
        url: 'pubs.cgi?year=' + year + '&type=' + type,
        dataType:'html',
        data: {action: 'request'},
        success: function(data){
            \$('#result').html(data);
        },
        eroor: function(){
            alert('Error');
        }
    });
}
</script>";


print $query->end_html;
exit(0);
