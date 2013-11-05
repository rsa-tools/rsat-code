----------------------------------------------------------------
-- given a list of gene identifiers, get similar genes (results of
-- BLAST comparisons) from the database

SELECT
	gi,
	database_organism,
	query_gi,
	query_organism,
	raw_score,
	length,
	bits,
	significance,
	rank
	
FROM
	blast_hits
WHERE
	query_gi = 16131776
