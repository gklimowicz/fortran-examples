FORTRAN_LANG_URL = https://fortran-lang.org/en/packages
FIND_IS_FORTRAN= -type f \( -iname "*.f" -o -iname "*.f[0-9][0-9]" -o -iname "*.ftn" \)

all:	all-files all-fortran-files all-projects all-fortran-files-lc stats


# Create a list of all Fortran projects we have,
# all files, and all Fortran files.
all-projects:	Makefile
	/bin/ls -d */ | sed -e 's;/$$;;' >"$@"

all-files:	. Makefile
	find * ! -path ".git*" -a -type f -print \
	| sort >"$@"

all-fortran-files:	. Makefile
	find * ${FIND_IS_FORTRAN} -print \
	| sort >"$@"

# Word count each Fortran file
all-fortran-files-lc:	all-fortran-files Makefile
	cat all-fortran-files | tr '\n' '\0' \
	| xargs -0 wc -l \
	| grep -v 'total$$' >"$@"

# Word count each Fortran file
all-projects-lc:	all-fortran-files Makefile
	for P in `cat all-projects`; do \
		echo "$$P\t$$(find $$P ${FIND_IS_FORTRAN} -print0 \
			| xargs -0 cat \
			| wc -l)"; \
	done >"$@"

# Print some moderately interesting stats about the repositories.
# stats: all-files all-fortran-files all-projects-lc Makefile
stats: FORCE
	@printf "%'11d projects\n" $$(wc -l <all-projects)
	@printf "%'11d files\n" $$(wc -l <all-files)
	@printf "%'11d Fortran files\n" $$(wc -l <all-fortran-files)
	@printf "%'11d Fortran lines\n" $$(awk '{ sum += $$1} END { print sum}' <all-fortran-files-lc)
	@printf "%'11d cpp directive lines\n" $$(cat all-fortran-files | tr '\n' '\0' | xargs -0 cat | grep '^#' | wc -l)
	@echo "Largest projects:"; \
	 sort -rn -k 2 all-projects-lc | head -10 | awk $$'{ printf "  %-10s %\'10d\\n", $$1, $$2}'


# Add new projects we may find lying about.
add-new-projects:	fortran-lang-new-projects

update-existing-projects:
	git submodule foreach git pull

# Get ne projects for fortran-lang.org.
# Only look at Git projects (net sourceforge or whatever).
fortran-lang-new-projects:	fortran-lang-projects
	for P in `cat fortran-lang-projects`; do \
		case "$$P" in git*) ;; *) continue;; esac; \
		if [[ ! -d "$$(basename "$$P")" ]]; then \
			: echo "$$P exists already"; \
		else \
			echo "$$P new"; \
			git submodule add ssh://git@$$P.git \
			|| git submodule add https:$$P.git; \
		fi; \
	done

fortran-lang-projects:	fortran-lang-category-urls
	for P in `cat fortran-lang-category-urls`; do \
		curl -L --no-progress-meter "$$P" \
		| sed -n -e '/^<p class="rubric"><span class="/s;.*href="https*://\([^"]*\).*;\1;p'; \
	done \
	| sort -u >"$@"

fortran-lang-category-urls:	FORCE
	curl -L --no-progress-meter "${FORTRAN_LANG_URL}" 2>&1 \
	| sed -n -e '/^<h2><a class="reference internal" href="\([^"]*\).*/s;;${FORTRAN_LANG_URL}/\1;p' \
	| sort -u >"$@"

FORCE:
