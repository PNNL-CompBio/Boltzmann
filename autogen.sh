#!/bin/sh

aclocal &&
automake --copy --force-missing --add-missing &&
autoconf --force

