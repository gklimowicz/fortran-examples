.SUFFIXES:
.PHONY: update help

#
# modelE directory structure
MODEL_E_ROOT = .
MODEL_DIR = $(MODEL_E_ROOT)/model
CONFIG_DIR = $(MODEL_E_ROOT)/config

ifeq ($(VERBOSE_OUTPUT),NO)
  MAKEFLAGS=-s
endif

CURRENT_RELEASE=$(shell [ -s $(MODEL_DIR)/version ] && cat $(MODEL_DIR)/version)

default:
	@echo "All make commands should be executed from 'model/decks' "
	@echo "directory. Please 'cd decks' before running any commands"

help:
	$(MAKE) -C $(MODEL_DIR) help

sinclude $(CONFIG_DIR)/rules.mk

update:
	@if [ "$(RELEASE)" = "" ] ; then \
	  echo "You have to specify new release with RELEASE=..." ; exit 1 ; \
	fi
	@echo "--- Starting automatic update from CVS repository ---"
	@echo "Checking if current directory tree is a branch"
	@if cvs status | grep "Sticky *Tag:.*none" ; then \
	  echo "Some of your files are on the main trunk. Will not update."; \
	  echo "You need to be on a branch to do automatic updates." ;\
	  echo "Exiting" ; exit 1 ; \
	fi
	@echo "Checking if the model knows its current release"
	@if [ "$(CURRENT_RELEASE)" = "" ] ; then \
	  echo "No information on current release. Exiting." ; exit 1 ; \
	fi
	@if cvs status -v README | grep $(CURRENT_RELEASE) ; then \
	  echo "Current release: $(CURRENT_RELEASE)" ; \
	else \
	  echo "No tag $(CURRENT_RELEASE) in repository." ;\
	  echo "Will not update. Exiting"; exit 1 ; \
	fi
	@echo "Searching repository for specified tag: $(RELEASE)"
	@if cvs status -v README | grep $(RELEASE) ; then \
	  echo "Found tag: $(RELEASE)" ; \
	else \
	  echo "The tag $(RELEASE) is not present in CVS repository" ;\
	  echo "Will not update. Exiting"; exit 1 ; \
	fi
	@echo "Checking current directory tree"
	@if cvs status | grep "Status: *Needs" ; then \
	  echo "The files listed above are not up-to-date" ; \
	  echo "Do 'cvs update' first. Exiting." ; exit 1; \
	fi
	@echo "Commiting your latest changes to CVS"
	@cvs commit -m "last commit before update to $(RELEASE)" ; \
	if [ ! $$? = 0 ] ; then \
	  echo "Something went wrong with commit. Aborting..." ; exit 1 ; \
	fi
	@echo "Just in case, last check before we start the update."
	@if cvs status | grep "Status:" | grep -v "Up-to-date" ; then \
	  echo "The files listed above are not Up-to-date" ; \
	  echo "Can not do automatic update. Aborting." ; exit 1 ; \
	fi
	@echo "----------   updating   ----------"
	@echo "(using: -j $(CURRENT_RELEASE) -j $(RELEASE) )"
	@cvs update -j $(CURRENT_RELEASE) -j $(RELEASE) ; \
	if [ ! $$? = 0 ] ; then \
	  echo "Something went wrong with update..." ; \
	  echo "You have to check it manually. Sorry about that ..." ; \
	  exit 1 ; \
	else \
	  echo "--- Successfully updated to $(RELEASE) ---" ; \
	fi
	@echo "Checking for possible conflicts"
	@if cvs status | grep "Status:.*conflict" ; then \
	  echo "The files listed above had conflicts on merge." ; \
	  echo "Please check them. " ; \
	fi
	@echo "After you check your model for possible conflicts"
	@echo "It will be a good idea to commit this new version to cvs"
	@echo "with 'cvs commit' "
	@echo "----------   finished   ----------"





















