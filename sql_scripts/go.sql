----------------------------------------------------------------
-- select relationships between terms, filtered on child term name
select
--	term1_id,
--	term2_id,
	t1.id parent_id,
	t1.name parent_name,
	t2.id child_id,
	t2.name child_name,
	relationship_type_id
from 
	term2term,
	term t1,
	term t2 
where
--	t1.id = 7374
	t2.name like 'methionine metabolism%' 
--and
--	relationship_type_id=5
and 
	term1_id=t1.id 
and 
	term2_id=t2.id
limit 100
;
