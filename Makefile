##############################################################################
#   PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
#   Copyright (C) 2007-2009 University of Tennessee / UT-Battelle
#
#   Code Author:  Doug Hyatt
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

CC=	gcc

SRC=	main.c \
	gene.c \
	dprog.c \
	node.c \
	sequence.c \
	training.c \
	bitmap.c

PROC=	prodigal

CFLAGS=	-O3 -Wall

LIBS=	-lm
LDFLAGS=	$(LDIRS) $(LIBS)


OBJS=	${SRC:.c=.o}


all:	$(PROC)

clean:	$(OBJS)
	/bin/rm -rf $(OBJS) $(PROC)

$(PROC):	$(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS) $(CFLAGS)

bitmap.o:	bitmap.c
	$(CC) -c $(CFLAGS) bitmap.c
dprog.o:	dprog.c
	$(CC) -c $(CFLAGS) dprog.c
gene.o:	gene.c
	$(CC) -c $(CFLAGS) gene.c
main.o:	main.c
	$(CC) -c $(CFLAGS) main.c
node.o:	node.c
	$(CC) -c $(CFLAGS) node.c
sequence.o:	sequence.c
	$(CC) -c $(CFLAGS) sequence.c
training.o:	training.c
	$(CC) -c $(CFLAGS) training.c
