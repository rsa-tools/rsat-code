################################################################
#
# This script creates users for the mysql database
# It assumes to be logged as mysql su


################################################################
# Create databases for amaze parsing result
create database swissprot;
GRANT ALL ON swissprot.* TO swissprot@localhost;
GRANT ALL ON swissprot.* TO jvanheld@localhost;
GRANT ALL ON swissprot.* TO rekins@localhost;

create database kegg;
GRANT ALL ON kegg.* TO kegg@localhost IDENTIFIED BY kegg;
GRANT ALL ON kegg.* TO jvanheld@localhost;
GRANT ALL ON kegg.* TO rekins@localhost;

create database sigtrans;
GRANT ALL ON sigtrans TO sigtrans@localhost;
GRANT ALL ON sigtrans TO jvanheld@localhost;
GRANT ALL ON sigtrans TO rekins@localhost;


################################################################
# For info: commands required to throw a user
################################################################

use mysql;

DELETE FROM user WHERE user='kegg';
DELETE FROM tables_priv WHERE user='kegg';
DELETE FROM db WHERE user='kegg';

DELETE FROM user WHERE user='swissprot';
DELETE FROM tables_priv WHERE user='swissprot';
DELETE FROM db WHERE user='swissprot';

DELETE FROM user WHERE user='sigtrans';
DELETE FROM tables_priv WHERE user='sigtrans';
DELETE FROM db WHERE user='sigtrans';

DELETE FROM user WHERE user='jvanheld';
DELETE FROM tables_priv WHERE user='jvanheld';
DELETE FROM db WHERE user='jvanheld';

DELETE FROM user WHERE user='rekins';
DELETE FROM tables_priv WHERE user='rekins';
DELETE FROM db WHERE user='rekins';

drop database swissprot;
drop database kegg;
drop database sigtrans;
