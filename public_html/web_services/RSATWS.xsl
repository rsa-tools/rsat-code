<?xml version='1.0' encoding='UTF-8'?>
<xsl:stylesheet version='1.0'
xmlns:soap='http://schemas.xmlsoap.org/wsdl/soap/'
xmlns:wsdl='http://schemas.xmlsoap.org/wsdl/'
xmlns:xsl='http://www.w3.org/1999/XSL/Transform'
xmlns:fo='http://www.w3.org/1999/XSL/Format' 
xmlns:xsd='http://www.w3.org/2001/XMLSchema'
xmlns:html='http://www.w3.org/1999/xhtml'>

  <xsl:template match='/'>
    <xsl:apply-templates/>
  </xsl:template>
  
  <xsl:template match='wsdl:definitions'>
    <html>
      <head>
        <title>Documentation generated from WSDL</title>
      </head>	
      <body>
        <table border='0' align='left' width='80%'>
          <tr><th>
              <table border='1' align='left' width='100%' cellpadding='10'>
                <xsl:apply-templates select='wsdl:portType'/>
                <xsl:apply-templates select='wsdl:service'/>
                <xsl:apply-templates select='wsdl:portType/wsdl:operation'/>
              </table>
          </th></tr>
        </table>
      </body>
    </html>		
  </xsl:template>

  <xsl:template match='wsdl:portType'>
    <tr><th bgcolor='#000080'>
        <font color='#FFFFFF' size='5'><xsl:value-of select='@name'/></font>
    </th></tr>
  </xsl:template>

  <xsl:template match='wsdl:service'>
    <tr><th align='left'>
        <table width='100%' cellpadding='10'>
          <tr>
            <th width='200' align='left' valign='top'><font color='#2200FF' size='4'>Service Documentation</font></th>
            <th><xsl:value-of select='wsdl:documentation'/></th>
          </tr>
          <tr>
            <th align='left' valign='top'><font color='#2200FF' size='4'>Server Address</font></th>
            <th><xsl:value-of select='wsdl:port/soap:address/@location'/></th>
          </tr>
        </table>
    </th></tr>
  </xsl:template>

  <!-- Template invoked for each method -->
  <xsl:template match='wsdl:portType/wsdl:operation'>
    <xsl:variable name = 'varInputName' select = 'wsdl:input/@name' />
    <xsl:variable name = 'varOutputName' select = 'wsdl:output/@name' />    
        
    <tr><th align='left' valign='top'>
        <table width='100%' cellpadding='5' border='0' rules='none' frame='void' >
          <tr align='left' valign='top'>
            <th width='200' align='left' valign='top'><font size='4' color='#888888'>Method</font></th>
            <td align='left' valign='top'><font size='4' color='#222288'><xsl:value-of select='@name'/></font></td>
          </tr>
          <tr align='left' valign='top'>
            <th align='left' valign='top'><font color='#888888' size='3'>Description</font></th>
            <td align='left' valign='top'><font color='#444444' size='3'><xsl:value-of select='wsdl:documentation'/><br/></font></td>
          </tr>
          <tr align='left' valign='top'>
            <th align='left' valign='top'><font color='#888888' size='3'>Parameters</font></th>
            <td align='left' valign='top'><font color='#444444' size='3'>
                <font color='#0000DD'>Input Parameters</font>
                <xsl:apply-templates select='ancestor::*//wsdl:types/xsd:schema/xsd:complexType[@name=$varInputName]/xsd:sequence/xsd:element'/><br/>
                <font color='#0000DD'>Output Parameters</font>
                <xsl:apply-templates select='ancestor::*//wsdl:types/xsd:schema/xsd:complexType[@name=$varOutputName]/xsd:sequence/xsd:element'/><br/></font>
            </td>
          </tr>
        </table>	
    </th></tr>
  </xsl:template>

  <!-- Template invoked for each parameter of a method -->	
  <xsl:template match='xsd:element'>
    <tr align='left'>
      <td height='35'>
        <font color='#228888'><xsl:value-of select='@name'/></font><br/>
	<xsl:if test='@minOccurs >= 1'>(required)</xsl:if>
      </td>
      <td>
	<pre><!-- font face='Times' size='+1' --><span style="font-family:Times; font-size:12pt;color:black"><xsl:value-of select='xsd:annotation/xsd:documentation'/><!-- /font --></span></pre><!-- br/ -->
	<i>
	  <xsl:if test='@type = "xsd:string"'>Type = String</xsl:if>
	  <xsl:if test='@type = "xsd:int"'>Type = Integer</xsl:if>
      	  <xsl:if test='@type = "xsd:boolean"'>Type = Boolean</xsl:if>
      	  <xsl:if test='@type = "xsd:float"'>Type = Float</xsl:if>
	</i>
      </td>
    </tr>
  </xsl:template>

</xsl:stylesheet>
