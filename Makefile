SHELL = /bin/bash

# We compress the largest files to keep GitHub
# from complaining about file sizes.
BUILD_FILES=all-projects.txt \
	    all-files.txt all-files.txt.gz \
	    all-fortran-files.txt \
	    all-fortran-files-attr.txt all-fortran-files-attr.txt.gz \
	    all-fortran-files-type.txt all-fortran-files-type.txt.gz \
	    all-fortran-files-lc.txt \
	    all-fortran-files-fixed.txt \
	    all-fortran-files-free.txt \
	    all-projects-lc.txt \
	    all-projects-fortran-file-count.txt \
	    stats.txt

# Update the .gz file only if it changes.
%.txt.gz: %.txt
	gzip --keep "$<"

all: ${BUILD_FILES}

# Update projects from their remote repos.
update:
	git submodule update --remote

# Create a list of all Fortran projects we have,
# all files, and all Fortran files.
all-projects.txt: origins.txt .gitmodules
	(git submodule foreach -q 'echo $$sm_path'; \
	 sed -n -e 's/#.*//' -e 's/:.*//p' <origins.txt) \
    | sort >"$@"

all-files.txt: all-projects.txt
	find `cat all-projects.txt` ! -path ".git*" -a -type f -print \
	| sort >"$@"

all-files.txt.gz: all-files.txt

# Some files that match the fortran-file-patterns.txt list are not Fortran files.
# Eliminate them explicitly.
all-fortran-files.txt:	all-files.txt fortran-file-patterns.txt \
			fortran-exceptions.txt
	grep -i -f fortran-file-patterns.txt all-files.txt \
	| fgrep -v -f fortran-exceptions.txt >"$@"
	wc -l "$@"

# Compute the attributes in parallel and concatenate results.
# This is awkward, but cuts the time to generate them down dramatically.
# Note that since all-fortran-files.txt is sorted, the attributes
# file will be as well, without having to explicitly sort it.
# We touch it at the end, for the unusual circumstance
# where the file looks out of date, but no changes were
# made via `cpif,.
all-fortran-files-attr.txt:	all-fortran-files.txt bin/determine-attributes
	bin/create-stats-file all-fortran-files.txt bin/determine-attributes \
	| bin/cpif "$@"
	wc -l "$@"

all-fortran-files-attr.txt.gz: all-fortran-files-attr.txt

# All fixed-form Fortran files
all-fortran-files-fixed.txt:	all-fortran-files-attr.txt
	gawk  -F '\t' '/form:fixed/ { print $$1 }' \
	      <"$<" >"$@"
	wc -l "$@"

# All free-form Fortran files
all-fortran-files-free.txt:	all-fortran-files-attr.txt
	gawk -F '\t' '/form:free/ { print $$1 }' \
	      <"$<" >"$@"
	wc -l "$@"

# Results of running "file" on each Fortran file.
# This is useful for finding non-Fortran files that are
# disguised with Fortran-file-like names.
all-fortran-files-type.txt:	all-fortran-files.txt /usr/bin/file
	bin/create-stats-file \
		all-fortran-files.txt /usr/bin/file \
	| bin/cpif "$@"
	touch "$@"
	wc -l "$@"

all-fortran-files-type.txt.gz: all-fortran-files-type.txt

# Line count of Fortran files in each project
all-fortran-files-lc.txt:	all-fortran-files-attr.txt
	gawk -F '\t' \
	    '{ file = $$1; \
	         lc = gensub(/.*lines:([0-9]*).*/, "\\1", 1, $$0)+0; \
	         printf("%8d\t%s\n", lc, file); \
	     }' \
	      <all-fortran-files-attr.txt >"$@"
	wc -l "$@"

all-projects-lc.txt:   all-fortran-files-lc.txt
	gawk -F '\t' \
	    'function new_proj() { \
	         printf("%8d %s\n", lc, last_proj); \
	         last_proj = project; \
	         lc = 0; \
	     } \
	     NR == 1 { \
	         last_proj = gensub(/\/.*/, "", 1, $$2); \
	         lc = $$1; \
	         next; \
	     } \
	     { \
	         project = gensub(/\/.*/, "", 1, $$2); \
	         if (project != last_proj) \
	             new_proj(); \
	         lc += $$1; \
	     } \
	     END { \
	         new_proj(); \
	     }' \
	      <all-fortran-files-lc.txt >"$@"
	wc -l "$@"

# List number of Fortran files in each project
all-projects-fortran-file-count.txt:	all-fortran-files-attr.txt
	gawk -F '\t' \
	     'BEGIN { last_proj = "" } \
	      { proj = gensub(/\/.*/, "", 1, $$1); \
	        if (proj != last_proj) { \
	             if (last_proj != "") \
	                 printf "%s\t%d\n", last_proj, sum; \
	             last_proj = proj; \
	             sum = 1; \
	        } else { \
	             sum += 1; \
	        } \
	      } \
	      END { printf "%s\t%d\n", last_proj, sum }' \
	      <all-fortran-files-attr.txt >"$@"
	wc -l "$@"

# Print some moderately interesting stats about the repositories.
stats.txt:  all-projects.txt all-projects-lc.txt all-files.txt \
			all-fortran-files.txt all-fortran-files-lc.txt \
			all-projects-fortran-file-count.txt
	@(printf "%'12d projects\n" $$(wc -l <all-projects.txt); \
	  printf "%'12d files\n" $$(wc -l <all-files.txt); \
	  printf "%'12d Fortran files\n" $$(wc -l <all-fortran-files.txt); \
	  printf "%'12d Fortran lines\n" \
		$$(gawk -F '\t' \
	            '{ sum += gensub(/.*lines:([0-9]*).*/, "\\1", 1, $$0)+0 } \
                    END { print sum }' \
	            <all-fortran-files-attr.txt); \
	  printf "%'12d Fortran 77 lines\n" \
		$$(gawk -F '\t' \
	            '/form:fixed/ { \
                         sum += gensub(/.*lines:([0-9]*).*/, "\\1", 1, $$0)+0} \
                     END { print sum }' \
		<all-fortran-files-attr.txt); \
	  printf "%'12d cpp directive lines\n" \
	        $$(gawk -F '\t' \
	            '/cpreprocessor/ { \
                         sum += gensub(/.*cpreprocessor:([0-9]*).*/, "\\1", 1, $$0)+0 \
	             } \
                     END { print sum }' \
	           <all-fortran-files-attr.txt); \
	  echo "Largest projects:"; \
	  sort -rn -k 1 all-projects-lc.txt \
		| head -10 \
		| awk $$'{ printf "%\'12d %s\\n", $$1, $$2}'; \
	  echo "Largest Fortran files:"; \
	  sort -rn -k 1 all-fortran-files-lc.txt \
		| head -10 \
		| awk $$'{ printf "%\'12d %s\\n", $$1, $$2}') \
	| tee "$@"
	wc -l "$@"


# Add new projects we may find lying about.
add-new-projects:	fortran-lang-new-projects beliavsky-new-projects

# Get new projects for fortran-lang.org.
# Only look at Git projects (not sourceforge or whatever).
fortran-lang-new-projects:	fortran-lang-projects.txt
	for P in `cat fortran-lang-projects.txt`; do \
		case "$$P" in \
			git*) ;; \
			*) continue;; \
		esac; \
		D="$$(echo "$$P" | awk -F / '{ print $$3 "@" $$2 }')"; \
	    if grep "^$$D:" project-exceptions.txt; then \
			echo "$$D is on the project exception list."; \
			echo ""; \
		elif [[ -d $$D ]]; then \
			echo "$$D exists already"; \
		else \
			if git submodule add "ssh://git@$$P" "$$D"; then \
			    echo "$$D new"; \
			    git submodule update --init "$$D"; \
		    else \
		        echo "$$D failed; probably no longer available"; \
		    fi; \
		fi; \
	done

fortran-lang-projects.txt:	fortran-lang-category-urls.txt
	for P in `cat fortran-lang-category-urls.txt`; do \
		curl -L --no-progress-meter "$$P" \
		| sed -n -e '/^<h2 class="rubric".*<span class="/s;.*href="https*://\([^"]*\).*;\1;p'; \
	done \
	| sort -u >"$@"


# Get Fortran projects from https://fortran-lang.org/
# and add new ones as submodules of this project.
# These are organized into separate pages based on
# project category, so we gather the category pages first.

FORTRAN_LANG_URL = https://fortran-lang.org/en/packages
fortran-lang-category-urls.txt:	FORCE
	curl -L --no-progress-meter "${FORTRAN_LANG_URL}" 2>&1 \
	| sed -n -e '/^<h2><a class="reference internal" href="\([^"]*\).*/s;;${FORTRAN_LANG_URL}/\1;p' \
	| sort -u >"$@"

# Get new projects from github.com/Beliavsky/Fortran-code-on-GitHub.
# Only look at Git projects (not sourceforge or whatever).
beliavsky-new-projects:	beliavsky-projects.txt
	for P in `cat beliavsky-projects.txt`; do \
	    case "$$P" in \
	        git*) ;; \
	        *) continue;; \
	    esac; \
	    D="$$(echo "$$P" | awk -F / '{ print $$3 "@" $$2 }')"; \
	    if grep "^$$D:" project-exceptions.txt; then \
	        echo "$$D is on the project exception list."; \
		echo ""; \
	    elif [[ -d $$D ]]; then \
	        echo "$$D exists already"; \
	    elif git submodule add "ssh://git@$$P" "$$D"; then \
	           echo "$$D new"; \
	           git submodule update --init "$$D"; \
	    else \
	        echo "$$D failed; probably no longer available"; \
	    fi; \
	done

BELIAVSKY_URL = https://github.com/Beliavsky/Fortran-code-on-GitHub
beliavsky-projects.txt:	FORCE
	curl -L --no-progress-meter "${BELIAVSKY_URL}" \
		| sed -n -e 's;^<p dir="auto"><a href="https://\([^"]*\)".*;\1;p' \
	        | sed -e 's/#.*//' \
		| sort -u >"$@"

clean: FORCE
	rm -f ${BUILD_FILES}

FORCE:
