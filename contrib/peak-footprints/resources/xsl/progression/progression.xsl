<?xml version="1.0" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="html" encoding="ISO-8859-1" />

<xsl:template match="/">
<html>
	 <xsl:apply-templates select="pipeline"/>
</html>
</xsl:template>


<xsl:template match="pipeline">
    <head>
        <title>
            progression of <xsl:value-of select="@name"/>
        </title>
    </head>

    <body style="font-family:Arial; font-size:12pt;">

	<div style="text-align:center; font-weight:bold; font-size:20pt; background-color:#6666cc; color:white; padding:4px"> 
		Workflow <xsl:value-of select="@name"/>
	</div>
        <br></br>
	<div style="text-align:center; font-weight:bold;color:black; padding:4px">
		Start time : <xsl:value-of select="@startTime"/>
	</div>
	<div style="text-align:center; font-weight:bold;color:black; padding:4px">
		<xsl:if test="@endTime != '0'">
			End time : <xsl:value-of select="@endTime"/>	
		</xsl:if>
	</div>
        <br></br>
        <br></br>
        <table ALIGN="center" WIDTH='80%'>
            <tr HEIGHT='5' BGCOLOR="black"> 
                <th> 
                    <div style="text-align:center; font-weight:bold;color:white; padding:4px"> 
                        Component Branch and Name
                    </div>
                </th>
                <th> 
                    <div style="text-align:center; font-weight:bold;color:white; padding:4px"> 
                        Elapsed Time
                    </div>
                </th>
                <th> 
                    <div style="text-align:center; font-weight:bold;color:white; padding:4px"> 
                        Status
                    </div>
                </th>
            </tr>
            <tr HEIGHT='5' FGCOLOR="white" BGCOLOR="white"> 
		<td> <div style="text-align:center; font-weight:bold;color:white; padding:4px">.</div> </td>
		<td> <div style="text-align:center; font-weight:bold;color:white; padding:4px">.</div> </td>
                <td> <div style="text-align:center; font-weight:bold;color:white; padding:4px">.</div> </td>
            </tr>
            
            <xsl:apply-templates select="component"/>
        </table>
    </body>
	
</xsl:template>



<xsl:template match="component">

    <xsl:choose>
        <xsl:when test="@status='Not Started'">
            <tr HEIGHT='5' BGCOLOR="grey">
                <td WIDTH='60%'>
                    <div style="margin-left:30px; font-weight:bold;color:white; padding:4px"> 
                        <xsl:value-of select="@rank"/> - <xsl:value-of select="@branch"/> <xsl:value-of select="@displayName"/>
                    </div>
                </td>
                <td WIDTH='20%'> 
                </td>
                <td WIDTH='20%'> 
                    <div style="text-align:center; color:white; padding:4px"> 
                         <xsl:value-of select="@status"/>
                    </div>
                </td>
            </tr>
            <tr HEIGHT='10' >
                <td WIDTH='60%'>
                </td>
                <td>
                </td>
            </tr>
        </xsl:when>
        
        <xsl:when test="@status='Running'">
            <tr HEIGHT='20' BGCOLOR="#6666cc">
                <td WIDTH='60%'>
                    <div style="margin-left:30px; font-weight:bold; color:white; padding:4px">
                        <xsl:value-of select="@rank"/> - <xsl:value-of select="@branch"/> <xsl:value-of select="@displayName"/>
                    </div>
                </td>
                
                <td WIDTH='20%'> 
                    <div style="text-align:center; color:white; padding:4px"> 
                        <xsl:value-of select="@elapsedTime"/>
                    </div>
                </td>
                
                <td WIDTH='20%'> 
                    <div style="text-align:center; color:white; padding:4px"> 
                        <xsl:value-of select="@status"/>
                    </div>
                </td>
            </tr>
	    <xsl:if test="@progression">
		    <tr HEIGHT='10'>
			<td WIDTH='60%' HEIGHT='30' >
			    <div style="margin-left:80px; font-weight:bold; color:black; padding:4px">
				Progression : <xsl:value-of select="@progression"/>
			    </div>
			    <br></br>
			</td>
			<td WIDTH='20%' HEIGHT='30' > 
			    <div style="color:black; padding:4px"> 

			    </div>
			</td>
			<td WIDTH='20%' HEIGHT='30' > 
			    <div style="color:black; padding:4px"> 

			    </div>
			</td>
		    </tr>
	   </xsl:if>

            <xsl:apply-templates select="task"/>
        </xsl:when>
        
        <xsl:when test="@status='Executed'">
            <tr HEIGHT='10' BGCOLOR="#55cc55">
                <td WIDTH='60%'>
                    <div style="margin-left:30px; font-weight:bold;color:white; padding:4px"> 
                        <xsl:value-of select="@rank"/> - <xsl:value-of select="@branch"/> <xsl:value-of select="@displayName"/>
                    </div>
                </td>
                <td WIDTH='20%'> 
                    <div style="text-align:center; color:white; padding:4px"> 
                        <xsl:value-of select="@elapsedTime"/>
                    </div>
                </td>
                <td WIDTH='20%'> 
                    <div style="text-align:center; color:white; padding:4px"> 
                        <xsl:value-of select="@status"/>
                    </div>
                </td>
            </tr>
            <tr HEIGHT='10' >
                <td WIDTH='60%'>
                </td>
                <td>
                </td>
            </tr>
            <xsl:apply-templates/>
        </xsl:when>

        <xsl:when test="@status='Resumed'">
            <tr HEIGHT='10' BGCOLOR="#55cc55">
                <td WIDTH='60%'>
                    <div style="margin-left:30px; font-weight:bold;color:white; padding:4px"> 
                        <xsl:value-of select="@rank"/> - <xsl:value-of select="@branch"/> <xsl:value-of select="@displayName"/>
                    </div>
                </td>
                <td WIDTH='20%'> 
                </td>
                <td WIDTH='20%'> 
                    <div style="text-align:center; color:white; padding:4px"> 
                        <xsl:value-of select="@status"/> 
                    </div>
                </td>
            </tr>
            <tr HEIGHT='10' >
                <td WIDTH='60%'>
                </td>
                <td>
                </td>
            </tr>
            <xsl:apply-templates/>
        </xsl:when>    
    
        <xsl:when test="@status='Failed'">
            <tr HEIGHT='10'  BGCOLOR="#cc5555">
                <td WIDTH='60%'>
                    <div style="margin-left:30px; font-weight:bold; color:white; padding:4px">
                        <xsl:value-of select="@rank"/> - <xsl:value-of select="@branch"/> <xsl:value-of select="@displayName"/>
                    </div>
                </td>
                <td WIDTH='20%'> 
                </td>
                <td WIDTH='20%'> 
                    <div style="text-align:center; color:white; padding:4px"> 
                        <xsl:value-of select="@status"/>
                    </div>
                </td>
            </tr>
            <tr HEIGHT='10' >
                <td WIDTH='60%'>
                </td>
                <td>
                </td>
            </tr>
        </xsl:when>
        <xsl:otherwise></xsl:otherwise>
    </xsl:choose>

</xsl:template>

<xsl:template match="task">

    <tr HEIGHT='10'>
        <td WIDTH='60%' HEIGHT='30' >
            <div style="margin-left:80px; font-weight:bold; color:black; padding:4px">
                <xsl:value-of select="@name"/>
		<xsl:if test="@progression = '0.0%'">
			...
		</xsl:if>
		<xsl:if test="@progression != '0.0%'">
			: <xsl:value-of select="@progression"/>
		</xsl:if>

            </div>
	    <br></br>
        </td>
        <td WIDTH='20%' HEIGHT='30' > 
            <div style="color:black; padding:4px"> 

            </div>
        </td>
        <td WIDTH='20%' HEIGHT='30' > 
            <div style="color:black; padding:4px"> 

            </div>
        </td>
    </tr>
    
</xsl:template>


</xsl:stylesheet>
