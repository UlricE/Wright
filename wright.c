#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <sdb.h>

static char *url, *table, *id, *birth;
static char *dam, *sire, *coefficient;
static char *kinship;
static char *kincondition = "true";
static int verbose = 0;
static int max_id;

struct animal {
	int id, sire, dam, sx, dx;
	int a1, a2;	/* alleles for founder alleles surviving simulation */
	int dmin, dmax;	/* min and max generation distance from founders */
	double ic;
};

static struct animal *animal;
static int animal_max;
static double **rc;
static int fas_count = 10000;
static int fas_total;
static int allele;
static int *alleles;

static int imax(int a, int b)
{
	if (a < b) return b;
	return a;
}

static int imin(int a, int b)
{
	if (a < b) return a;
	return b;
}

static void usage(void)
{
	printf("Usage:\n"
               "  wright [options]\n"
               "\n"
               "Options:\n"
               "  -u url         SDB url to access the database\n"
               "  -t table       The name of the table\n"
               "  -i id          The name of the primary key field\n"
               "  -b birth       The name of the birth field\n"
               "  -d dam         The name of the dam field\n"
               "  -s sire        The name of the sire field\n"
               "  -c coefficient The name of the coefficient field\n"
               "  -k kinship     The name of the mean kinship field\n"
               "  -n condition   Additional kinship selection argument\n"
               "  -f condition   Additional FAS selection argument\n"
               "  -l #           FAS simulation loop count\n"
               "  -m             SQL query to select which parts of kinship matrix to save\n"
	       "  -g             SQL query to pick out founder representation to save\n"
	       "  -D             Dump generational depth for each animal\n"
               "  -v             Verbose\n"
               "  -h             Display this information\n");
	printf("\n");
	printf("Example:\n");
	printf("  wright -u mysql:uid=myusername:pwd=mypassword:db=mydatabase \\\n"
		"\t-t register -i id -b birth -d dam_id -s sire_id \\\n"
		"\t-c ic \\\n"
		"\t-n \"g2008 != '' and bogus is null\" -k mk2008 \\\n"
		"\t-n \"g2009 != '' and bogus is null\" -k mk2009\n"
		"\t-m \"select r1.id, r2.id from kanin_g_register r1 join kanin_g_register r2 where r1.id < r2.id and r1.kon != r2.kon and r1.g2009 != '' and r2.g2009 != ''\"\n"
		"\t-g \"select r.id, f.id from kanin_g_register r join kanin_g_register f where r.g2009 != '' and f.mor_id is null\"\n"
		"\t-l 10000\n"
		"\t-f \"g2009 != '' and bogus is null\"\n");
	printf("\n");
	printf("Parameter order matters.\n");

	exit(0);
}

static int lookup_id(int id)
{
	int i;
	if (id == -1) return max_id;	/* founder */
	for (i = 0; i < max_id; i++) {
		if (animal[i].id == id) break;
	}
	return i;
}

static int collect_cb(int n, char **p, void *closure)
{
	if (max_id >= animal_max) {
		animal_max += 1000;
		animal = realloc(animal, animal_max*sizeof *animal);
		if (animal == NULL) {
			fprintf(stderr, "Out of memory allocating animal table\n");
			exit(EXIT_FAILURE);
		}
	}
	if (p[0] == NULL) {
		if (verbose) fprintf(stderr, "null id\n");
		return 1;	/* shouldn't happen */
	}
	animal[max_id].id = atoi(p[0]);
	/* use -1 to flag founders, i.e. sire/dam unknown */
	if (p[1] == NULL) {
		if (verbose) fprintf(stderr, "null sire\n");
		animal[max_id].sire = -1;
	} else {
		animal[max_id].sire = atoi(p[1]);
	}
	if (p[2] == NULL) {
		if (verbose) fprintf(stderr, "null dam\n");
		animal[max_id].dam = -1;
	} else {
		animal[max_id].dam = atoi(p[2]);
	}
	max_id++;
	return 0;
}

static int null_cb(int n, char **p, void *closure)
{
	return 0;
}

static double reco(int i1, int i2)
{
	if (i1 == max_id || i2 == max_id) return 0;

	return rc[i1][i2];
}

static void relationships(void)
{
	char query[1024];
	int i, j, n;
	static int done = 0;

	if (verbose) fprintf(stderr, "Calculating relationships\n");
	if (done) {
		if (verbose) fprintf(stderr, "Already done\n");
		return;
	}

	if (table == NULL) {
		fprintf(stderr, "Table name missing\n");
		exit(EXIT_FAILURE);
	}
	if (id == NULL) {
		fprintf(stderr, "Name of id field missing\n");
		exit(EXIT_FAILURE);
	}
	if (birth == NULL) {
		fprintf(stderr, "Name of birth field missing\n");
		exit(EXIT_FAILURE);
	}
	if (dam == NULL) {
		fprintf(stderr, "Name of dam field missing\n");
		exit(EXIT_FAILURE);
	}
	if (sire == NULL) {
		fprintf(stderr, "Name of sire field missing\n");
		exit(EXIT_FAILURE);
	}

	snprintf(query, sizeof query,
		"select %s, %s, %s\n"
		"from %s\n"
		"order by %s\n",
		id, sire, dam, table, birth);
	if (verbose) fprintf(stderr, "SDB: '%s'\n", query);
	n = sdb_query(url, query, collect_cb, NULL);

	if (verbose) fprintf(stderr, "sdb_query returns %d\n", n);

	rc = malloc(max_id * sizeof *rc);
	if (rc == NULL) {
		fprintf(stderr, "Out of memory allocating relationship table\n");
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < max_id; i++) {
		animal[i].sx = lookup_id(animal[i].sire);
		animal[i].dx = lookup_id(animal[i].dam);
		if (verbose >= 3) {
			fprintf(stderr, "%d: <%d,%d,%d,%d,%d>\n", i, animal[i].id,
				animal[i].sire, animal[i].dam,
				animal[i].sx, animal[i].dx);
		}
		rc[i] = malloc(max_id * sizeof **rc);
		if (rc[i] == NULL) {
			fprintf(stderr, "Out of memory allocating relationship table\n");
			exit(EXIT_FAILURE);
		}
	}

	if (verbose) fprintf(stderr, "Fill out the entire relationship matrix\n");
	for (i = 0; i < max_id; i++) {
		/* an animal's relationship coefficient to itself is */
		/* 1 (because it is completely related to itself */
		/* plus half the relationship coefficient of its parents */
		rc[i][i] = 1+reco(animal[i].sx, animal[i].dx)/2;
		/* the relationship coefficient of two different animals is */
		/* half the rc of the first animal and the other's sire */
		/* plus half the rc of the first animal and the other's dam */
		for (j = i+1; j < max_id; j++) {
			rc[i][j] = (reco(i, animal[j].sx) +
				    reco(i, animal[j].dx))/2;
			/* of course, j relates to i as i relates to j */
			rc[j][i] = rc[i][j];
		}
		if (verbose >= 3) fprintf(stderr, "%d: %f\n", i, (rc[i][i]-1)*100);
	}
	if (verbose) fprintf(stderr, "Assigning unique founder alleles\n");
	allele = 1;
	for (i = 0; i < max_id; i++) {
		/* founders do not have parents */
		if (animal[i].sx == max_id && animal[i].dx == max_id) {
			animal[i].a1 = allele++;
			animal[i].a2 = allele++;
		}
	}
	alleles = malloc(allele*sizeof *alleles);
	if (alleles == NULL) {
		fprintf(stderr, "Out of memory allocating alleles table\n");
		exit(EXIT_FAILURE);
	}
	if (verbose) fprintf(stderr, "Assigning min and max generation depth\n");
	for (i = 0; i < max_id; i++) {
		/* founders have depth = 0 */
		if (animal[i].sx == max_id && animal[i].dx == max_id) {
			animal[i].dmin = animal[i].dmax = 0;
		} else {
			int sx = animal[i].sx, dx = animal[i].dx;
			animal[i].dmin = 1+imin(animal[sx].dmin, animal[dx].dmin);
			animal[i].dmax = 1+imax(animal[sx].dmax, animal[dx].dmax);
		}
	}
#if 0	/* can't do this now */
	for (i = 0; i < max_id; i++) free(rc[i]);
	free(rc);
#endif
	done = 1;
}

static void inbreeding(void)
{
	int i;
	char query[1024];

	if (verbose) fprintf(stderr, "Calculating inbreeding; results in '%s'\n", coefficient);
	relationships();

	for (i = 0; i < max_id; i++) {
		if (fabs(rc[i][i]) < 0.000001) {
			snprintf(query, sizeof query,
				"update %s\n"
				"set %s = 0\n"
				"where %s = %d",
				table, coefficient, id, animal[i].id);
		} else {
			snprintf(query, sizeof query,
				"update %s\n"
				"set %s = %f\n"
				"where %s = %d",
				table, coefficient, rc[i][i], id, animal[i].id);
		}
		if (verbose >= 2) fprintf(stderr, "%s\n", query);
		sdb_query(url, query, null_cb, NULL);
	}
}

static struct kinship {
	int id;
	int index;
	double mk;
} *kinship_table;

static int max_kinship;

static int kinship_cb(int n, char **p, void *closure)
{
	if (verbose >= 2) fprintf(stderr, "kinship_cb(%d, '%s', %p\n", n, p[0], closure);
	kinship_table[max_kinship].id = atoi(p[0]);
	/* the index is into the rc table */
	kinship_table[max_kinship].index = lookup_id(kinship_table[max_kinship].id);
	if (verbose >= 2) fprintf(stderr, "%d: index=%d <=> id=%d\n", max_kinship,
				kinship_table[max_kinship].index,
				kinship_table[max_kinship].id);
	max_kinship++;
	return 0;
}

static int fas_cb(int n, char **p, void *closure)
{
	int id1 = atoi(p[0]);
        int i1 = lookup_id(id1);
	alleles[animal[i1].a1]++;
	alleles[animal[i1].a2]++;
	return 0;
}

static void fas_simulation(char *fas_condition)
{
	char query[1024];
	int i, j, acount;
	if (verbose) fprintf(stderr, "Simulating Founder Alleles Surviving for '%s', %d loops\n", fas_condition, fas_count);
	relationships();
	fas_total = 0;
	for (j = 0; j < fas_count; j++) {
		for (i = 0; i < allele; i++) {
			alleles[i] = 0;
		}
		for (i = 0; i < max_id; i++) {
			if (animal[i].sx != max_id) {	/* not founder */
				if (random() & 1) {
					animal[i].a1 = animal[animal[i].sx].a1;
				} else {
					animal[i].a1 = animal[animal[i].sx].a2;
				}
				if (random() & 1) {
					animal[i].a2 = animal[animal[i].dx].a1;
				} else {
					animal[i].a2 = animal[animal[i].dx].a2;
				}
			}
		}
		snprintf(query, sizeof query,
			"select id from %s where %s",
			table, fas_condition);
		if (verbose >= 2) fprintf(stderr, "%s\n", query);
		sdb_query(url, query, fas_cb, NULL);
		acount = 0;
		for (i = 0; i < allele; i++) {
			if (alleles[i]) acount++;
		}
		if (verbose) fprintf(stderr, "Simulation %d, allele count: %d\n", j, acount);
		fas_total += acount;
	}
	printf("Total Alleles Surviving average of %d simulations: %.02f\n", fas_count, (double)fas_total/fas_count);
}

static void depth(void)
{
	int i;

	if (verbose) fprintf(stderr, "Dumping generation distance from founders\n");
	relationships();
	for (i = 0; i < max_id; i++) {
		printf("Depth:\t%d\t%d\t%d\n", animal[i].id, animal[i].dmin, animal[i].dmax);
	}
}

static void mean_kinship(void)
{
	char query[1024];
	int i, j;

	if (verbose) fprintf(stderr, "Calculating mean kinship where %s, results in '%s'\n", kincondition, kinship);

	relationships();

	kinship_table = malloc(max_id*sizeof *kinship_table);
	if (kinship_table == NULL) {
		fprintf(stderr, "Out of memory allocating kinship table\n");
		exit(EXIT_FAILURE);
	}
	max_kinship = 0;
	snprintf(query, sizeof query,
		"select %s\n"
		"from %s\n"
		"where %s\n",
		id, table, kincondition);
	if (verbose) fprintf(stderr, "%s\n", query);
	sdb_query(url, query, kinship_cb, NULL);

	snprintf(query, sizeof query,
		"update %s set %s = null\n",
		table, kinship);
	if (verbose >= 2) fprintf(stderr, "%s\n", query);
	sdb_query(url, query, null_cb, NULL);

	for (i = 0; i < max_kinship; i++) {
		int index_i = kinship_table[i].index;
		double total_kinship = 0;
		for (j = 0; j < max_kinship; j++) {
			int index_j = kinship_table[j].index;
#if 0	/* apparently the definition includes kinship with self */
			if (i != j) total_kinship += rc[index_i][index_j];
#else
			total_kinship += rc[index_i][index_j];
#endif
		}
		/* I do not know why the 2 is necessary */
		kinship_table[i].mk = total_kinship/max_kinship/2;
		snprintf(query, sizeof query,
			"update %s\n"
			"set %s = %f\n"
			"where %s = %d",
			table, kinship, kinship_table[i].mk, id, kinship_table[i].id);
		if (verbose >= 2) fprintf(stderr, "%s\n", query);
		sdb_query(url, query, null_cb, NULL);
	}
	free(kinship_table);
}

static int matrix_cb(int n, char **p, void *closure)
{
        int id1 = atoi(p[0]), id2 = atoi(p[1]);
        int i1, i2;

        if (id1 == 0 || id2 == 0) return 1;     // can't do that
        i1 = lookup_id(id1);
        i2 = lookup_id(id2);
        if (i1 == max_id || i2 == max_id) return 1;    // not that either
        printf("Kinship: %d %d %f\n", id1, id2, rc[i1][i2]);
	return 0;
}

static int founders_cb(int n, char **p, void *closure)
{
        int id1 = atoi(p[0]), id2 = atoi(p[1]);
        int i1, i2;

        if (id1 == 0 || id2 == 0) return 1;     // can't do that
        i1 = lookup_id(id1);
        i2 = lookup_id(id2);
        if (i1 == max_id || i2 == max_id) return 1;    // not that either
        if (rc[i1][i2] > 0) printf("Founders: %d %d %f\n", id1, id2, rc[i1][i2]);
	return 0;
}

static void save_kinship(char *query)
{
	int n;
        if (verbose) fprintf(stderr, "%s\n", query);
        n = sdb_query(url, query, matrix_cb, NULL);
        if (verbose) fprintf(stderr, "sdb_query returns %d\n", n);
}

static void save_founders(char *query)
{
	int n;
	if (verbose) fprintf(stderr, "%s\n", query);
	n = sdb_query(url, query, founders_cb, NULL);
	if (verbose) fprintf(stderr, "sdb_query returns %d\n", n);
}

int main(int argc, char **argv)
{
	int c;
	char *opt = "u:t:i:b:d:s:c:k:n:m:g:l:f:vhD";

	while ((c = getopt(argc, argv, opt)) != -1) {
		switch (c) {
		case 'u':
			url = sdb_open(optarg);
			break;
		case 't':
			table = optarg;
			break;
		case 'i':
			id = optarg;
			break;
		case 'b':
			birth = optarg;
			break;
		case 'd':
			dam = optarg;
			break;
		case 's':
			sire = optarg;
			break;
		case 'c':
			coefficient = optarg;
			inbreeding();	/* update the database */
			break;
		case 'k':
			kinship = optarg;
			mean_kinship();	/* update the database */
			break;
		case 'n':
			kincondition = optarg;
			break;
		case 'm':
			save_kinship(optarg);
			break;
		case 'g':
			save_founders(optarg);
			break;
		case 'l':
			fas_count = atoi(optarg);
			break;
		case 'f':
			fas_simulation(optarg);
			break;
		case 'v':
			verbose++;
			break;
		case 'D':
			depth();
			break;
		default:
			usage();
		}
	}

	if (url) sdb_close(url);

	return 0;
}
