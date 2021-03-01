<?xml version="1.0" encoding="utf-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text"/>
<xsl:template match="/">

</xsl:template>

<xsl:template match="/">
<xsl:apply-templates select="root"/>
</xsl:template>


<xsl:template match="root">
<xsl:text>entropy	exp_Mcount	exp_rpkm	exp_total	full_rpkm	gene	idis_metadata	is_sample	project_desc	sample_id	source_name	sra_id	taxid	var	
</xsl:text>
<xsl:apply-templates select="doc"/>
</xsl:template>

<xsl:template match="doc">
<xsl:value-of select="field[@name='entropy']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='exp_Mcount']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='exp_rpkm']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='exp_total']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='full_rpkm']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='gene']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='id']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='is_metadata']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='is_sample']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='project_desc']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='sample_id']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='source_name']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='sra_id']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='taxid']"/>
<xsl:text>	</xsl:text>
<xsl:value-of select="field[@name='var']"/>
<xsl:text>
</xsl:text>
</xsl:template>


</xsl:stylesheet>
