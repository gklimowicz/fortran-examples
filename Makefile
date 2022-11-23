FORTRAN_LANG_URL = https://fortran-lang.org/en/packages

all:	all-fortran-files all-projects all-fortran-files-wc stats

# Find all Fortran files in all these repositories.
all-fortran-files:	. Makefile
	find * -type f \( -iname "*.f" -o -iname "*.f[0-9][0-9]" -o -iname "*.ftn" \) -print \
	| sort >"$@"

# Create a list of all Fortran projects we have here.
all-projects:	Makefile
	/bin/ls -d */ | sed -e 's;/$$;;' >"$@"

# Word count each Fortran file
all-fortran-files-wc: all-fortran-files Makefile
	cat all-fortran-files | tr '\n' '\0' \
	| xargs -0 wc -l \
	| grep -v 'total$$' >"$@"

# Print some moderately interesting stats about the repositories.
stats: all-fortran-files Makefile
	@echo "$$(cat all-fortran-files | tr '\n' '\0' | xargs -0 cat | wc -l | tr -d ' ') total lines"
	@echo "$$(cat all-fortran-files | tr '\n' '\0' | xargs -0 cat | grep '^#' | wc -l | tr -d ' ') directive lines"


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

fortran-lang-projects: fortran-lang-category-urls
	for P in `cat fortran-lang-category-urls`; do \
		curl -L --no-progress-meter "$$P" \
		| sed -n -e '/^<p class="rubric"><span class="/s;.*href="https*://\([^"]*\).*;\1;p'; \
	done \
	| sort -u >"$@"

fortran-lang-category-urls: FORCE
	curl -L --no-progress-meter "${FORTRAN_LANG_URL}" 2>&1 \
	| sed -n -e '/^<h2><a class="reference internal" href="\([^"]*\).*/s;;${FORTRAN_LANG_URL}/\1;p' \
	| sort -u >"$@"

FORCE:
