#!/bin/bash

REV = $(shell git rev-parse --verify HEAD)
CC = g++

a.out: main.h main.cpp
	@rm -rf revision.h
	@echo "string revision = \"$(REV)\";" > revision.h
	$(CC) -O3 -g -ffast-math main.cpp
	@rm -rf revision.h

clean:
	rm -rf *.txt
	rm -rf *.log
	rm -rf a.out
