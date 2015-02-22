
PACKAGE = wright
VERSION = 0.0.3beta1

PREFIX = /usr/local
BINDIR = $(PREFIX)/bin

SOURCES = wright.c
DOCS = README INSTALL ChangeLog COPYING AUTHORS TODO

CFLAGS = -O -Wall -g

wright : wright.c
	$(CC) `sdb-config --cflags --libs` wright.c -o wright `sdb-config --cflags --libs`

install: wright
	cp wright $(BINDIR)

uninstall:
	rm -f $(BINDIR)/wright

dist:
	rm -rf $(PACKAGE)-$(VERSION)
	mkdir $(PACKAGE)-$(VERSION)
	cp $(SOURCES) Makefile $(DOCS) $(PACKAGE)-$(VERSION)
	tar cf - $(PACKAGE)-$(VERSION) | gzip > $(PACKAGE)-$(VERSION).tar.gz
	rm -rf $(PACKAGE)-$(VERSION)

clean:
	$(RM) wright *.o *.bak *~

