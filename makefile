#!/bin/bash

REV = $(shell git rev-parse --verify HEAD)

exe:
	@rm -rf revision.h
	@echo "string revision = \"$(REV)\";" > revision.h
	g++ -O3 main.cpp
	@rm -rf revision.h
