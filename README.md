# Wright

This program calculates Wright's coefficient of inbreeding, using
data from a database.

It also calculates mean kinship for each animal in a population.

Example:

wright -u mysql:uid=myusername:pwd=mypassword:db=mydatabase \
	-t register -i id -b birth -d dam_id -s sire_id \
	-c ic \
	-n "g2008 != ''" -k mk2008 \
	-n "g2009 != ''" -k mk2009

This will first calculate the relationship matrix for all animals in the
database. Then inbreeding coefficients will be stored for all animals.
Finally two sets of mean kinships will be calculated and stored. Any number
of calculations can be performed. Being able to do so is useful because
the values change every time an animal is added or removed.

Parameter order is significant.
