<?xml version="1.0" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
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
        <title>
            Results of <xsl:value-of select="@name"/>
        </title>
    </head>

    <body style="font-family:Arial; font-size:12pt;">

	<div style="text-align:center; font-size:20pt; background-color:#555555; color:white; padding:4px"> 
		<b>Results of <xsl:value-of select="@name"/></b>
	</div>


        <br></br>
        <br></br>

	<!-- ................................................................................................. -->
	<!-- Table containg the information on data input and global statistics                                -->
	<!-- ................................................................................................. -->
	<table ALIGN="center">
		<tr><td><b>Reference Species</b></td><td>:<xsl:value-of select="@referenceSpecies"/></td></tr>
		<tr><td><b>Reference Motif</b></td><td>:<xsl:value-of select="@referenceMotif"/></td></tr>
		<tr><td><b>Aligned Species</b></td><td>:<xsl:value-of select="@alignedSpecies"/></td></tr>
		<tr><td><b>Number of initial sequences</b></td><td>:<xsl:value-of select="@bedSequenceNumber"/></td></tr>
		<tr><td><b>Number of sequences associated to MSA</b></td><td>:<xsl:value-of select="@bedSequenceNumberWithMSA"/></td></tr>
		<tr><td><b>Number of detected conserved regions</b></td><td>:<xsl:value-of select="@conservedRegionNumber"/></td></tr>
	</table>


	<br></br>
        <br></br>


	<!-- ................................................................................................. -->
	<!-- Table containg the information on sequences, conserved regions and conserved blocks distributions -->
	<!-- ................................................................................................. -->
	<div style="text-align:center; font-size:10pt; color:black; padding:4px"> 
	<table ALIGN="center" WIDTH='50%'>
		<tr>
			<td ALIGN="center" >
				<a>
					<xsl:attribute name="href"> <xsl:value-of select="@BEDSequencesSizeHistogramGraph"/> </xsl:attribute>
					<img>
						<xsl:attribute name="src"> <xsl:value-of select="@BEDSequencesSizeHistogramGraph"/> </xsl:attribute>
						<xsl:attribute name="title"> Distribution of sequences size</xsl:attribute>
						<xsl:attribute name="border"> 0 </xsl:attribute>
						<xsl:attribute name="height"> 120 </xsl:attribute>
					</img>
				</a>
			</td>
			<td ALIGN="center" >
			 	<a> 
					<xsl:attribute name="href"> <xsl:value-of select="@MSASizeHistogramGraph"/> </xsl:attribute>
					<img>
						<xsl:attribute name="src"> <xsl:value-of select="@MSASizeHistogramGraph"/> </xsl:attribute>
						<xsl:attribute name="title"> Distribution of MSA size</xsl:attribute>
						<xsl:attribute name="border"> 0 </xsl:attribute>
						<xsl:attribute name="height"> 120 </xsl:attribute>
					</img>
				</a>
			</td>
			<td ALIGN="center" >
				<a> 
					<xsl:attribute name="href"> <xsl:value-of select="@ConservedRegionsSizeGraph"/> </xsl:attribute>
					<img>
						<xsl:attribute name="src"> <xsl:value-of select="@ConservedRegionsSizeGraph"/> </xsl:attribute>
						<xsl:attribute name="title"> Distribution of conserved regions size</xsl:attribute>
						<xsl:attribute name="border"> 0 </xsl:attribute>
						<xsl:attribute name="height"> 120 </xsl:attribute>
					</img>
				</a>
			</td>

		</tr>
		<tr HEIGHT='5' BGCOLOR="#ffffff">
			<th WIDTH='33%'>
				Distribution of sequences sizes
			</th>
			<th WIDTH='33%'>
				Size distribution of conserved regions under peaks
			</th>
			<th WIDTH='33%'> 
				Size distribution of conserved blocks under peaks
			</th>

		</tr>
	</table>
	</div>



	<!-- ................................................................................................. -->
	<!-- Table headers containg a line for each motif and its statistics data                              -->
	<!-- ................................................................................................. -->

        <br></br>
        <br></br>

	<div style="text-align:center; font-size:10pt; color:black; padding:4px"> 
        <table ALIGN="center" WIDTH='95%' border="1">
            <tr HEIGHT='5' BGCOLOR="#339933"> 
                <th><b>Name</b></th>
                <th><b>Hits number</b></th>
                <th><b>Expected hits number</b></th>
                <th><b>Ratio with hits numbers</b></th>
                <th><b>Hypergeometric P-value</b></th>
                <th><b>Chi2 P-value</b></th>
		<th><b>Chi2</b></th>
                <th>% Overlapping sites with <xsl:value-of select="@referenceMotif"/></th>
                <th ALIGN="left"><b> 
                        <ul>
				<li style='margin-left:20px'>Family</li>
				<li style='margin-left:20px'>Class</li>
				<li style='margin-left:20px'>Type</li>
			</ul>
                </b></th>
                <th><b>Logo</b></th>
                <th><b>Motif hit distribution in peaks</b></th>
                <th><b>Distance to nearest <xsl:value-of select="@referenceMotif"/></b></th>
                <!--th><b>Peak scores of hits</b></th-->
           </tr>
            
	    <!-- Apply the template at the family level -->
            <xsl:apply-templates select="family/motif"/>

        </table>
	</div>


	<br></br>
	<br></br>



	<!-- ................................................................................................. -->
	<!-- Documentation                                                                                     -->
	<!-- ................................................................................................. -->

	<center></center>
	<table ALIGN="center" WIDTH='95%'>
		<tr ALIGN="center" BGCOLOR="#339933"><td colspan="2"><h2>Documentation</h2></td></tr>
		<tr ALIGN="center" HEIGHT='5' BGCOLOR="#66cc66">
			<th>
				<b>Table entry</b>
			</th>
			<th>
				<b>Details</b>
			</th>
		</tr>
		<tr><td WIDTH='30%'>Reference Species</td><td>The species on which the initial data where analyzed from</td></tr>
		<tr><td WIDTH='30%'>Reference Motif</td><td>The motif that have been used for the ChIP-seq analysis</td></tr>
		<tr><td WIDTH='30%'>Aligned Species</td><td>The list of species that vae been used to create the Multiple Sequence Alignements (MSA)</td></tr>
		<tr><td WIDTH='30%'>Number of initial sequences</td><td>The initial number of sequence found in the input data file</td></tr>
		<tr><td WIDTH='30%'>Number of sequences associated to MSA</td><td>The number of sequence that were effectively associate to a MSA</td></tr>
		<tr><td WIDTH='30%'>Number of detected conserved regions</td><td>The number of conserved regions identified through all the MSA</td></tr>
	</table>
	<br></br>
	<table ALIGN="center" WIDTH='95%'>
		<tr><td WIDTH='30%'>Distribution of sequences sizes </td><td>The histogram of the sizes of the initial sequences</td></tr>
		<tr><td WIDTH='30%'>Distribution of MSA sizes </td><td>The histogram of the sizes of the composed MSA</td></tr>
		<tr><td WIDTH='30%'>Distribution of conserved regions sizes </td><td>The histogram of the sizes of the identified conserved regions</td></tr>
	</table>
	<br></br>
	<table ALIGN="center" WIDTH='95%'>
		<tr><td WIDTH='30%'>Name</td><td>The name of the family or of the motif</td></tr>
		<tr><td WIDTH='30%'>Hits Number</td><td>The number of occurence of the motif through the conserved regions</td></tr>
		<tr><td WIDTH='30%'>Expected hits number/Ratio with hits numbers</td><td>The number of occurence of the motif permuted matrix/The ratio of this number with the motif hits number</td></tr>
		<tr><td WIDTH='30%'>Hypergeometric P-value</td><td>The p-value computed from the hypergeometric test</td></tr>
		<tr><td WIDTH='30%'>Chi2 P-value/Chi2</td><td>The p-value computed from the Chi2 test executed against an homogen distribution/the Chi2 value of the test</td></tr>
		<tr><td WIDTH='30%'>% Overlapping sites with <xsl:value-of select="@referenceMotif"/></td><td>The percentage of motif site that overlap a site of the reference motif (at least 50% of overlap)</td></tr>
		<tr><td WIDTH='30%'>Class/Type</td><td>The motif class and type</td></tr>
		<tr><td WIDTH='30%'>Logo</td><td>The motif logo</td></tr>
		<tr><td WIDTH='30%'>Distance to peak maximum</td><td>Histogram of the distances between each motif site and the sequence maximum. This maximum is the sequence center or the sequence maximal signal position. The corresponding homogen distributiopn is also graphed.</td></tr>
		<tr><td WIDTH='30%'>Distance to nearest <xsl:value-of select="@referenceMotif"/></td><td>Histogram of the distances between each motif site and the nearest reference motif site. Distances above the provided limit are ignored</td></tr>
	</table>
    </body>
	
</xsl:template>


<!-- ************************************************************************************************* -->
<!-- TEMPLATE ON THE MOTIF LEVEL                                                                       -->
<!-- ************************************************************************************************* -->

<xsl:template match="family/motif">
	<xsl:sort select="rank" order="ascending"/>
    <tr HEIGHT='10' ALIGN="center">
        <td WIDTH='15%' HEIGHT='30' >
		   <a>
			<xsl:attribute name="href"> <xsl:value-of select="@matrix"/> </xsl:attribute>
			<xsl:value-of select="@id"/> <br/>
			(<xsl:value-of select="@name"/>)
		   </a>
        </td>
        <td WIDTH='5%' HEIGHT='30' >
		<xsl:value-of select="@hitscore"/>       
        </td>
        <td WIDTH='5%' HEIGHT='30' >
		<xsl:value-of select="@hyphitsscore"/>       
        </td>
        <td WIDTH='10%' HEIGHT='30'> 
		<xsl:variable name="A"><xsl:value-of select="@hitscore"/></xsl:variable>
		<xsl:variable name="B"><xsl:value-of select="@hyphitsscore"/></xsl:variable>
		<xsl:value-of select="round(($A div $B)*100)"/>%
        </td>
        <td WIDTH='15%' HEIGHT='30' > 
		<xsl:value-of select="@hyppvalue"/>     
        </td>
        <td WIDTH='10%' HEIGHT='30' > 
		<xsl:value-of select="@chi2pvalue"/>     
        </td>
        <td WIDTH='10%' HEIGHT='30' > 
		<xsl:value-of select="@chi2"/>     
        </td>
        <td WIDTH='10%' HEIGHT='30' > 
		<xsl:value-of select="@ratioHomoLocation"/>
        </td>
        <td WIDTH='10%' HEIGHT='30' ALIGN="left"> 
		<ul>
			<!--li style='margin-left:-20px'><font size="2"><xsl:value-of select="@family"/></font></li-->
			<li style='margin-left:-20px'><font size="2"><xsl:value-of select="@class"/></font></li>
			<li style='margin-left:-20px'><font size="2"><xsl:value-of select="@type"/></font></li>
		</ul>
        </td>
        <td WIDTH='25%' HEIGHT='30' > 
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
        <td WIDTH='25%' HEIGHT='30' > 
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
        <td WIDTH='25%' HEIGHT='30' >
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
    </tr>
    
</xsl:template>


</xsl:stylesheet>
