################################################################
#
# This script creates users for the mysql database
# It assumes to be logged as mysql su


################################################################
# Create databases for amaze parsing result
create database swissprot;
GRANT ALL ON swissprot TO jvanheld;
GRANT ALL ON swissprot TO rekins;

create database kegg;
GRANT ALL ON kegg TO jvanheld;
GRANT ALL ON kegg TO rekins;

create database sigtrans;
GRANT ALL ON sigtrans TO jvanheld;
GRANT ALL ON sigtrans TO rekins;


################################################################
# For info: commands required to throw a user
################################################################
# use mysql;
DELETE FROM user WHERE user='jvanheld';
DELETE FROM tables_priv WHERE user='jvanheld';
DELETE FROM db WHERE user='jvanheld';

DELETE FROM user WHERE user='rekins';
DELETE FROM tables_priv WHERE user='rekins';
DELETE FROM db WHERE user='rekins';

drop database swissprot;
drop database kegg;
drop database sigtrans;
