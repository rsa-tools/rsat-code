----------------------------------------------------------------
-- given a list of gene names, get the IDs and description
SELECT
	n.id,
	n.names,
	g.id,
	g.name,
	g.taxid,
	g.organism
FROM
	cds_names n,
	cds g
WHERE
	n.names = 'metJ'
    AND
	n.id = g.id
;
