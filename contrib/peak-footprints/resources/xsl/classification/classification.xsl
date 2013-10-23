<?xml version="1.0" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<!--xsl:import href="JS-SortTable.xsl"/-->
<xsl:output method="html" encoding="ISO-8859-1" />


<!-- ************************************************************************************************* -->
<!-- HTML PAGE                                                                                         -->
<!-- ************************************************************************************************* -->


<xsl:template match="/">
<html>
	 <xsl:apply-templates select="classification"/>
</html>
</xsl:template>


<!-- ************************************************************************************************* -->
<!-- TEMPLATE ON THE ROOT LEVEL                                                                        -->
<!-- ************************************************************************************************* -->
<xsl:template match="classification">
    <head>
    	<!-- Loading the Jquery library -->
		<script src="http://code.jquery.com/jquery-1.10.1.min.js"></script>
		<!-- Load the JQuery DataTables library -->
		<script type="text/javascript" language="javascript" src="./jquery.dataTables.js"></script>
		<script type="text/javascript" charset="utf-8">
		    //Define function to show and hide parameter table
			$( "#hide_optional_parameters" ).click(function() {
	  			$( "#optional_parameters" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#show_optional_parameters" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_optional_parameters" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
			});
			
			$( "#show_optional_parameters" ).click(function() {
	  			$( "#optional_parameters" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#show_optional_parameters" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_optional_parameters" ).show( "fast", function() {
	    			// Animation complete.
	  			});
			});
		
			// Insure the DOM Is loaded
			$(document).ready(function() {
			//	$( "#optional_parameters" ).hide( "fast", function() {
	    		//	// Animation complete.
	  		//	});
	  			$( "#show_optional_parameters" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_optional_parameters" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  		});
	  		
	  		//Define function to show and hide global statistics table
			$( "#hide_global_statistics" ).click(function() {
	  			$( "#global_statistics" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#show_global_statistics" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_global_statistics" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
			});
			
			$( "#show_global_statistics" ).click(function() {
	  			$( "#global_statistics" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#show_global_statistics" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_global_statistics" ).show( "fast", function() {
	    			// Animation complete.
	  			});
			});
		
			// Insure the DOM Is loaded
			$(document).ready(function() {
			//	$( "#global_statistics" ).hide( "fast", function() {
	    		//	// Animation complete.
	  		//	});
	  			$( "#show_global_statistics" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_global_statistics" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  		});
	  		
	  		//Define function to show and hide documentation table
			$( "#hide_documentation" ).click(function() {
	  			$( "#documentation" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#show_documentation" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_documentation" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
			});
			
			$( "#show_documentation" ).click(function() {
	  		//	$( "#documentation" ).show( "fast", function() {
	    		//	// Animation complete.
	  		//	});
	  			$( "#show_documentation" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_documentation" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
			});
		
			// Insure the DOM Is loaded
			$(document).ready(function() {
				$( "#documentation" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#show_documentation" ).show( "fast", function() {
	    			// Animation complete.
	  			});
	  			$( "#hide_documentation" ).hide( "fast", function() {
	    			// Animation complete.
	  			});
	  		});
	  		
		
			/* Define two custom functions (asc and desc) for scientific number sorting */
			jQuery.fn.dataTableExt.oSort['scientific-number-case-asc'] = function(a,b) {
		          var x = parseFloat(a);
		          var y = parseFloat(b);
			      return ((x &lt; y) ? -1 : ((x &gt; y) ?  1 : 0));
			};
 
			jQuery.fn.dataTableExt.oSort['scientific-number-case-desc'] = function(a,b) {
		          var x = parseFloat(a);
		          var y = parseFloat(b);
			      return ((x &lt; y) ?  1 : ((x &gt; y) ? -1 : 0));
			};
			
			$(document).ready(function() {
				// Apply the JQuery DataTables style on the 'result_table' table
				$('#result_table').dataTable({
					// Sorting table with the third column (rank)
			 		"aaSorting": [[ 2, "asc" ]],
					"sPaginationType": "full_numbers",
					"aoColumns": [
            					null,
            					null,
            					null,
            					null,
            					null,
            					null,
					        	{ "sType": 'scientific-number-case' },
					        	{ "sType": 'scientific-number-case' },
				                null,
				                null,
				                null,
				                null,
				                null,
				                null
        				],
					"fnRowCallback": function( nRow, aData, iDisplayIndex ) {
						// Colorize the table row of the reference motif
						var motifname = aData[0];
						var refmotifname = "(" + '<xsl:value-of select="@referenceMotif"/>' + ")";
						if ( motifname.indexOf( refmotifname) !=-1 )
						{
							nRow.className='referenceMotifRow';
						}
						return nRow;
					}
				});
			} );
		</script>
    	<link rel="stylesheet" type="text/css" href = "./peak-footprints.css" media="screen"/>
        <title>
            Results of <xsl:value-of select="@name"/>
        </title>
    </head>

    <body style="font-family:Arial; font-size:12pt;">

	<table class="title_table">
	<!-- div style="text-align:center; font-size:20pt; background-color:#A9E2F3; color:#2E9AFE; padding:4px"-->
		<th>  
			Results of <xsl:value-of select="@name"/>
		</th>
	</table>


        <br></br>
        <br></br>

	<!-- ................................................................................................. -->
	<!-- Table containing the information on data input                                                      -->
	<!-- ................................................................................................. -->
	<table class="section_table">
		<tr><th>Analysis parameters</th></tr>
	</table>
	<table class="no_line_table" id="experiment_parameters">
		<tr><td><b>Reference Species</b></td><td> : <xsl:value-of select="@referenceSpecies"/></td></tr>
		<tr><td><b>Reference Motif</b></td><td> : <xsl:value-of select="@ReferenceMotifID"/> (<xsl:value-of select="@referenceMotif"/>)</td></tr>
		<tr><td><b>Aligned Species</b></td><td> : <xsl:value-of select="@alignedSpecies"/></td></tr>
		<tr><td><b>Input BED file</b></td><td> : <a><xsl:attribute name="href"><xsl:value-of select="@BEDFile"/> </xsl:attribute>Input BED File</a></td></tr>
		<tr><td><b>Output BED file</b></td><td> : 
			<a><xsl:attribute name="href"> <xsl:value-of select="@bedOutput"/> </xsl:attribute>BED Output</a> (
			<a>
				<xsl:attribute name="href"> 
					http://genome.ucsc.edu/cgi-bin/hgTracks?db=<xsl:value-of select="@referenceSpecies"/>&amp;hgct_customText=track%20type=bigBed%20name=<xsl:value-of select="@name"/>%20description=<xsl:value-of select="@name"/>%22%20visibility=full%20itemRgb=On%20bigDataUrl=http://pf:4test@139.124.66.4/pf/results/data/<xsl:value-of select="@name"/>/<xsl:value-of select="@name"/>_9_FinalOutputProcessor/BigBEDOutput/<xsl:value-of select="@name"/>_Motifs.bb
				</xsl:attribute>
				<xsl:attribute name="target">
					_blank
				</xsl:attribute>
				Show in UCSC Genome Browser
			</a>)
		</td></tr>
	</table>
	     
	<table class="section_table" width="50%">
		<tr><th id="show_optional_parameters">+ Additional parameters</th></tr>
		<tr><th id="hide_optional_parameters">- Additional parameters</th></tr>
	</table>
	<table class="no_line_table" id="optional_parameters">
	    <tr><td><b>Residue Conservation Limit</b></td><td> : <xsl:value-of select="@ResiduConservationLimit"/></td></tr>
	    <tr><td><b>Detection Window size</b></td><td> : <xsl:value-of select="@WindowSize"/></td></tr>
	    <tr><td><b>Detection Window conservation limit</b></td><td> : <xsl:value-of select="@WindowConservationLimit"/></td></tr>
	    <tr><td><b>Motif Databases root path</b></td><td> : <xsl:value-of select="@MotifDatabasePath"/></td></tr>
	    <tr><td><b>Motif Databases List</b></td><td> : 
	    	<xsl:value-of select="@MotifDatabaseFileList"/> 
	    	(<xsl:value-of select="@MotifDatabaseFileListSize"/> motifs)</td></tr>
	    <tr><td><b>Custom Motif Databases List</b></td><td> :
		    <xsl:variable name="CustomDB"><xsl:value-of select="@CustomMotifDatabaseFile"/></xsl:variable> 
	    	<xsl:if test="not($CustomDB = '')">
	    		<a><xsl:attribute name="href"> <xsl:value-of select="@CustomMotifDatabaseFile"/></xsl:attribute>Custom Motifs</a>
	    	 	(<xsl:value-of select="@CustomMotifDatabaseFileSize"/> motifs)
	    	</xsl:if>
	    	<xsl:if test="$CustomDB = ''">
	    		No Custom Motif
	    	</xsl:if>
	    </td></tr>
	    <tr><td><b>Max Hypergeometric p-value</b></td><td> : <xsl:value-of select="@MaxHypergeometricEValue"/></td></tr>
	    <tr><td><b>Max Chi2 p-value</b></td><td> : <xsl:value-of select="@MaxChi2EValue"/></td></tr>
	    <tr><td><b>Max Reported Motif Number</b></td><td> : <xsl:value-of select="@MaxMotifNumber"/></td></tr>
    </table>

	<!-- ................................................................................................. -->
	<!-- Table containing the information on global statistics                                               -->
	<!-- ................................................................................................. -->
	
	<table class="section_table">
		<tr id="show_global_statistics"><th>+ Global Statistics</th></tr>
		<tr id="hide_global_statistics"><th>- Global Statistics</th></tr>
	</table>
	<table class="statistic_table" id="global_statistics">
		<tr>
			<th></th>
			<th>Number</th>
			<th>Min Size</th>
			<th>Max Size</th>
			<th>Mean Size</th>
			<th>Total Size</th>
			<th>Size distribution</th>
		</tr>
		<tr><td><b>Initial Sequences</b></td>	
			<td align="center"><xsl:value-of select="@bedSequencesNumber"/></td>
			<td align="center"><xsl:value-of select="@bedSequencesMinSize"/></td>
			<td align="center"><xsl:value-of select="@bedSequencesMaxSize"/></td>
			<td align="center"><xsl:value-of select="@bedSequencesMeanSize"/></td>
			<td align="center"><xsl:value-of select="@bedSequencesTotalSize"/></td>
			<td align="center" >
			    <center>
				<a>
					<xsl:attribute name="href"> <xsl:value-of select="@BEDSequencesSizeHistogramGraph"/> </xsl:attribute>
					<img>
						<xsl:attribute name="src"> <xsl:value-of select="@BEDSequencesSizeHistogramGraph"/> </xsl:attribute>
						<xsl:attribute name="title"> Distribution of sequences size</xsl:attribute>
						<xsl:attribute name="height"> 90 </xsl:attribute>
						<xsl:attribute name="width"> 120 </xsl:attribute>
					</img>
				</a>
			   </center>
			   <center>
				<a>
					<xsl:attribute name="href"> <xsl:value-of select="@BEDSequencesSizeHistogram"/> </xsl:attribute>
					<font size="2">	Histogram </font>
				</a>
			   </center>
			</td>
		</tr>
		<tr><td><b>Sequences associated to MSA</b></td>	
			<td align="center"><xsl:value-of select="@MSANumber"/></td>
			<td align="center"><xsl:value-of select="@MSAMinSize"/></td>
			<td align="center"><xsl:value-of select="@MSAMaxSize"/></td>
			<td align="center"><xsl:value-of select="@MSAMeanSize"/></td>
			<td align="center"><xsl:value-of select="@MSATotalSize"/></td>
			<td align="center">
			   <center>
				<a> 
					<xsl:attribute name="href"> <xsl:value-of select="@MSASizeHistogramGraph"/> </xsl:attribute>
					<img>
						<xsl:attribute name="src"> <xsl:value-of select="@MSASizeHistogramGraph"/> </xsl:attribute>
						<xsl:attribute name="title"> Distribution of MSA size</xsl:attribute>
						<xsl:attribute name="height"> 90 </xsl:attribute>
						<xsl:attribute name="width"> 120 </xsl:attribute>
					</img>
				</a>
			   </center>
			   <center>
				<a>
					<xsl:attribute name="href"> <xsl:value-of select="@MSASizeHistogram"/> </xsl:attribute>
					<font size="2">	Histogram </font>
				</a>
			   </center>
			</td>
		</tr>
		<tr><td><b>Detected conserved Blocks</b></td>	
			<td align="center"><xsl:value-of select="@conservedBlocksNumber"/></td>
			<td align="center"><xsl:value-of select="@conservedBlocksMinSize"/></td>
			<td align="center"><xsl:value-of select="@conservedBlocksMaxSize"/></td>
			<td align="center"><xsl:value-of select="@conservedBlocksMeanSize"/></td>
			<td align="center"><xsl:value-of select="@conservedBlocksTotalSize"/></td>
			<td align="center">
			    <center>
				<a> 
					<xsl:attribute name="href"> <xsl:value-of select="@ConservedBlocksSizeGraph"/> </xsl:attribute>
					<img>
						<xsl:attribute name="src"> <xsl:value-of select="@ConservedBlocksSizeGraph"/> </xsl:attribute>
						<xsl:attribute name="title"> Distribution of conserved regions size</xsl:attribute>
						<xsl:attribute name="height"> 90 </xsl:attribute>
						<xsl:attribute name="width"> 120 </xsl:attribute>

					</img>
				</a>
			   </center>
			   <center>
				<a>
					<xsl:attribute name="href"> <xsl:value-of select="@ConservedBlocksSizeHistogram"/> </xsl:attribute>
					<font size="2">	Histogram </font>
				</a>
			   </center>
			</td>
		</tr>
	</table>


	<!-- ................................................................................................. -->
	<!-- Table headers containing a line for each motif and its statistics data                              -->
	<!-- ................................................................................................. -->

	<table class="section_table">
		<tr><th>Detected Motif Classification</th></tr>
	</table>
	<br></br>
	
        <table id="result_table">
        	<thead>
            <tr> 
                <th>Name</th>
                <th>Family</th>
				<th>Rank</th>
                <th>Hits number</th>
                <th>Expected hits number</th>
                <th>Ratio with hits numbers</th>
                <th>Hypergeometric P-value</th>
                <th>Chi2 P-value</th>
				<th>Chi2 value</th>
                <th>% Overlapping sites with <xsl:value-of select="@referenceMotif"/></th>
                <th>Class/ Type</th>
                <th>Logo</th>
                <th>Distance to peak center</th>
                <th>Distance to nearest <xsl:value-of select="@referenceMotif"/></th>
           </tr>
           </thead>
           <tbody>
	    		<!-- Apply the template at the family motifs level -->
            	<xsl:apply-templates select="family/motif"/>
           </tbody>

        </table>
        
	<!-- ................................................................................................. -->
	<!-- Documentation                                                                                     -->
	<!-- ................................................................................................. -->

	<br></br>
	<br></br>
	<table class="section_table">
		<tr id="show_documentation"><th>+ Documentation</th></tr>
		<tr id="hide_documentation"><th>- Documentation</th></tr>
	</table>
	<table class="no_line_table" id="documentation">
		<tr><td>Reference Species</td><td> : The species on which the initial data where analyzed from</td></tr>
		<tr><td>Reference Motif</td><td> : The motif that have been used for the ChIP-seq analysis</td></tr>
		<tr><td>Aligned Species</td><td> : The list of species that vae been used to create the Multiple Sequence Alignements (MSA)</td></tr>
		<tr><td>Number of initial sequences</td><td> : The initial number of sequence found in the input data file</td></tr>
		<tr><td>Number of sequences associated to MSA</td><td> : The number of sequence that were effectively associate to a MSA</td></tr>
		<tr><td>Number of detected conserved regions</td><td> : The number of conserved regions identified through all the MSA</td></tr>
		<tr><td></td></tr>
		<tr><td>Distribution of sequences sizes </td><td> : The histogram of the sizes of the initial sequences</td></tr>
		<tr><td>Distribution of MSA sizes </td><td> : The histogram of the sizes of the composed MSA</td></tr>
		<tr><td>Distribution of conserved regions sizes </td><td> : The histogram of the sizes of the identified conserved regions</td></tr>
		<tr><td></td></tr>
		<tr><td>Name</td><td> : The name and ID of the motif</td></tr>
		<tr><td>Family</td><td> : The motif family name</td></tr>
		<tr><td>Rank</td><td> : The rank of the motif in the classification</td></tr>
		<tr><td>Hits Number</td><td> : The number of occurence of the motif through the conserved regions</td></tr>
		<tr><td>Expected hits number</td><td> : The number of occurence of the motif permuted matrix</td></tr>
		<tr><td>Ratio with hits numbers</td><td> : The ratio of this number with the motif hits number</td></tr>
		<tr><td>Hypergeometric P-value</td><td> : The p-value computed from the hypergeometric test</td></tr>
		<tr><td>Chi2 P-value</td><td> : The p-value computed from the Chi2 test executed against an homogen distribution</td></tr>
		<tr><td>Chi2</td><td> : The Chi2 value of the Chi2 test executed against an homogen distribution</td></tr>
		<tr><td>% Overlapping sites with <xsl:value-of select="@referenceMotif"/></td><td> : The percentage of motif site that overlap a site of the reference motif (at least 50% of overlap)</td></tr>
		<tr><td>Class/Type</td><td> : The motif class and type</td></tr>
		<tr><td>Logo</td><td> : The motif logo</td></tr>
		<tr><td>Distance to peak center</td><td> : Histogram of the distances between each motif site and the peak center. The corresponding homogen distributiopn is also graphed.</td></tr>
		<tr><td>Distance to nearest <xsl:value-of select="@referenceMotif"/></td><td> : Histogram of the distances between each motif site and the nearest reference motif site. Distances above the provided limit are ignored</td></tr>
	</table>
	<br></br>
	<br></br>
   </body>
	
</xsl:template>


<!-- ************************************************************************************************* -->
<!-- TEMPLATE ON THE FAMILY LEVEL                                                                      -->
<!-- ************************************************************************************************* -->

<xsl:template match="family">

	
	<tr>
		<td>
			<b>Motif Family : <xsl:value-of select="@name"/></b>
		</td>
	</tr>

	<!-- Apply the template at the motif level -->
	<xsl:apply-templates select="motif"/>
    

</xsl:template>



<!-- ************************************************************************************************* -->
<!-- TEMPLATE ON THE MOTIF LEVEL                                                                       -->
<!-- ************************************************************************************************* -->

<xsl:template match="motif">

    <tr>
        <td>
		   <a>
			<xsl:attribute name="href"> <xsl:value-of select="@matrix"/> </xsl:attribute>
			<xsl:value-of select="@id"/> <br/>
			(<xsl:value-of select="@name"/>)
		   </a>
        </td>
        <td>
        	<xsl:value-of select="@family"/>	
        </td>
        <td>
		<xsl:value-of select="@rank"/>       
        </td>
        <td>
		<xsl:value-of select="@hitscore"/>       
        </td>
        <td>
		<xsl:value-of select="@hyphitsscore"/>       
        </td>
        <td> 
		<xsl:variable name="A"><xsl:value-of select="@hitscore"/></xsl:variable>
		<xsl:variable name="B"><xsl:value-of select="@hyphitsscore"/></xsl:variable>
		<xsl:value-of select="round(($A div $B)*100)"/>%
        </td>
        <td > 
		<xsl:value-of select="@hyppvalue"/>     
        </td>
        <td > 
		<xsl:value-of select="@chi2pvalue"/>     
        </td>
        <td > 
		<xsl:value-of select="@chi2"/>     
        </td>
        <td > 
		<xsl:value-of select="@ratioHomoLocation"/>
        </td>
        <td> 
        <xsl:value-of select="@class"/>/ <xsl:value-of select="@type"/>
        </td>
        <td> 
		<xsl:choose >
		   <xsl:when test="@logo">
				
			   <center><a> 
				<xsl:attribute name="href"> <xsl:value-of select="@logo"/> </xsl:attribute>
				<img>
					<xsl:attribute name="src"> <xsl:value-of select="@logo"/> </xsl:attribute>
					<xsl:attribute name="title"> <xsl:value-of select="@id"/> </xsl:attribute>
					<xsl:attribute name="border"> 0 </xsl:attribute>
					<xsl:attribute name="height"> 90 </xsl:attribute>
				</img>
			   </a></center>
   		   </xsl:when>
   		   <xsl:otherwise>
			-
   		   </xsl:otherwise>
		</xsl:choose >
        </td>
        <td> 
		<xsl:choose >
		   <xsl:when test="@distanceHistogramGraph">
			   <center><a> 
				<xsl:attribute name="href"> <xsl:value-of select="@distanceHistogramGraph"/> </xsl:attribute>
				<img>
					<xsl:attribute name="src"> <xsl:value-of select="@distanceHistogramGraph"/> </xsl:attribute>
					<xsl:attribute name="title"> <xsl:value-of select="@id"/> </xsl:attribute>
					<xsl:attribute name="height"> 90 </xsl:attribute>
					<xsl:attribute name="width"> 120 </xsl:attribute>
				</img>
			   </a></center>
			   <center>
				<a>
					<xsl:attribute name="href"> <xsl:value-of select="@distanceHistogram"/> </xsl:attribute>
					<font size="2">	Histogram </font>
				</a>
			   </center>
   		   </xsl:when>
   		   <xsl:otherwise>
			-
   		   </xsl:otherwise>
		</xsl:choose >
        </td>
        <td>
		<xsl:choose >
		   <xsl:when test="@coLocationHistogramGraph">
			   <center><a> 
				<xsl:attribute name="href"> <xsl:value-of select="@coLocationHistogramGraph"/> </xsl:attribute>
				<img>
					<xsl:attribute name="src"> <xsl:value-of select="@coLocationHistogramGraph"/> </xsl:attribute>
					<xsl:attribute name="title"> <xsl:value-of select="@id"/> </xsl:attribute>
					<xsl:attribute name="height"> 90 </xsl:attribute>
					<xsl:attribute name="width"> 120 </xsl:attribute>
				</img>
			   </a></center>
			   <center>
				<a>
					<xsl:attribute name="href"> <xsl:value-of select="@coLocationHistogram"/> </xsl:attribute>
					<font size="2">	Histogram </font>
				</a>
			   </center>
   		   </xsl:when>
   		   <xsl:otherwise>
			-
   		   </xsl:otherwise>
		</xsl:choose >
        </td>
        <!--td WIDTH='25%' HEIGHT='30' > 
		<xsl:choose >
		   <xsl:when test="@peakScoreHistogramGraph">
			   <center><a> 
				<xsl:attribute name="href"> <xsl:value-of select="@peakScoreHistogramGraph"/> </xsl:attribute>
				<img>
					<xsl:attribute name="src"> <xsl:value-of select="@peakScoreHistogramGraph"/> </xsl:attribute>
					<xsl:attribute name="title"> <xsl:value-of select="@id"/> </xsl:attribute>
					<xsl:attribute name="height"> 90 </xsl:attribute>
				</img>
			   </a></center>
			   <center>
				<a>
					<xsl:attribute name="href"> <xsl:value-of select="@peakScoreHistogram"/> </xsl:attribute>
					<font size="2">	Histogram </font>
				</a>
			   </center>
   		   </xsl:when>
   		   <xsl:otherwise>
			-
   		   </xsl:otherwise>
		</xsl:choose >
        </td-->
    </tr>
    
</xsl:template>


</xsl:stylesheet>