#!/bin/sh
#
# Copyright (C) 2008-2015 Michele Martone
# 
# This file is part of librsb.
# 
# librsb is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# librsb is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with librsb; see the file COPYING.
# If not, see <http://www.gnu.org/licenses/>.

# This script should be run as a cron sheduled job.
# It is intended for the librsb developer usage.

# You should set the FROMADDR (From Address), TOADDR (To Address), and SMTPHOST (Smtp Host) environment variables ...
#FROMADDR=user@host
#TOADDR=user@host
#SMTPHOST=smtp.host.tld

FROMADDR="${FROMADDR:-nightly.sh-for-user}"
TOADDR="${TOADDR:-user@domain.tld}"
SMTPHOST="${SMTPHOST:-some-smtp.tld}"

CONFIGURE_OPTS="$@"
RSB_SVN_REPOSITORY="${RSB_SVN_REPOSITORY:-svn+ssh://user@host/repository/trunk}"
TMPDIR=/tmp/

MS="librsb build FROMADDRILURE report"
alias if_err='[ $? == 0 ] || '
alias mail_stuff="mutt -F /dev/null -a autogen.log -a env.log -a make.log -a config.log -s librsb-automated-build  $TOADDR -e 'set from=$FROMADDR;set sendmail=\"msmtp --from $FROMADDR --host=$SMTPHOST\"' < librsb.log"

alias check='if_err { touch config.log librsb.log ; svn info  >> librsb.log  ; mail_stuff ;  }'
alias fail="return -1"


einfo() { echo  -- "[!]" $@ ; }

info()  { echo  -- "[*]" $@ ; }

die()   { einfo $@ ; exit -1; } 

get_rsb()
{
	#true;
	svn --force export $RSB_SVN_REPOSITORY librsb 
	#mkdir -p librsb 
}

autogen_error()
{
	# here should go automatic reporting of the failing configure.ac ...
	einfo "please see the autogen.log file"
}

configure_error()
{
	# here should go automatic reporting of the failing configure/config.log ...
	einfo "please see the config.log file"
}

make_error()
{
	# here should go automatic reporting of the failing make ...
	einfo "please see the $MAKELOG file"
}

date_ymd() { date +%Y%m%d ; }

build_rsb()
{
	MAKELOG="make.`date_ymd`.log"
	AUTOGENLOG="autogen.`date_ymd`.log"
	LOG="librsb.log"
	touch $MAKELOG || fail
	touch $AUTOGENLOG || fail
	rm -f $LOG || fail
	touch $LOG || fail
	date >>  $LOG || fail
	ln $MAKELOG make.log 
	ln $AUTOGENLOG autogen.log 
	check
	sh autogen.sh 2>&1 | tee $AUTOGENLOG  
	check
	#|| autogen_error "error generating initial librsb scripts"
	./configure "${CONFIGURE_OPTS}"
	check
	#|| configure_error "error configuring librsb"
	make clean 2>&1 | tee $MAKELOG
	check
	#|| make_error "error in making clean librsb"
	make       2>&1 | tee $MAKELOG
	check
	#alias mail_stuff
	mail_stuff
	#|| make_error "error making librsb"
}

[ -z "$RSB_SVN_REPOSITORY" ] && die "no librsb repository specified ?"
[ -z "$TMPDIR" ] && die "no temporary directory specified ?"


true
PROGRAMS="msmtp mutt svn"
which ${PROGRAMS}
if_err die "error looking for programs (need all in: $PROGRAMS) "
cd $TMPDIR 
if_err die "error stepping in $TMPDIR"
get_rsb
if_err die "error getting sources"
cd librsb
if_err die "error stepping in librsb directory"
env > env.log 
if_err die "error getting environment"
build_rsb 
if_err die "error building rsb"
cd -
if_err die "error stepping out of librsb directory"
rm -fR librsb 
if_err die "error cleaning up librsb directory"

