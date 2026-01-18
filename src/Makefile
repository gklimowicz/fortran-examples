SHELL = /bin/bash

PROJECTS_DIR=../projects
BUILD_ROOT=../build

# If you want to work on a subset of the projects,
# add 'LETTERS="[a-f]"' when invoking 'make'.
LETTERS=

# We compress the largest files to keep GitHub
# from complaining about file sizes.
BUILD_FILES=${BUILD_ROOT}/all-projects.txt \
	    ${BUILD_ROOT}/all-files.txt \
	    ${BUILD_ROOT}/all-files.txt.gz \
	    ${BUILD_ROOT}/all-fortran-files.txt \
	    ${BUILD_ROOT}/all-fortran-files-attr.txt \
	    ${BUILD_ROOT}/all-fortran-files-attr.txt.gz \
	    ${BUILD_ROOT}/all-fortran-files-type.txt \
	    ${BUILD_ROOT}/all-fortran-files-type.txt.gz \
	    ${BUILD_ROOT}/all-fortran-files-lc.txt \
	    ${BUILD_ROOT}/all-fortran-files-fixed.txt \
	    ${BUILD_ROOT}/all-fortran-files-free.txt \
	    ${BUILD_ROOT}/all-projects-lc.txt \
	    ${BUILD_ROOT}/all-projects-fortran-file-count.txt \
	    ${BUILD_ROOT}/stats.txt

# Update the .gz file only if it changes.
%.txt.gz: %.txt
	gzip --keep --force --best "$<"

all: check-no-empty-dirs ${BUILD_FILES}

# Update projects from their remote repos.
update: check-no-empty-dirs
	for D in ${PROJECTS_DIR}/${LETTERS}*@*; do \
	     (echo "Updating $$D..." \
	      && (builtin cd "$$D" \
	          && git fetch --tags --prune --prune-tags --force) \
	      && git submodule update --remote "$$D"); \
	done

# Perform sanity checks on the corpus.
# If there are directories that have no entries,
# that's an indicator that the submodules did not
# populate correctly.
check-no-empty-dirs:
	@builtin cd "${PROJECTS_DIR}"; \
	for D in *@*/; do \
	    echo "$$D"/[a-zA-Z0-9_]*; \
	done \
	    | sed -n -e 's/\/.*[*].*//p' \
	    | grep '.' \
	&& exit 1 \
	|| exit 0

# Create a list of all Fortran projects we have,
# all files, and all Fortran files.
${BUILD_ROOT}/all-projects.txt: origins.txt ../.gitmodules
	(git submodule foreach -q 'echo $$sm_path' \
	    | sed -e "s;$$(basename "${PROJECTS_DIR}")/;;"; \
	 sed -n -e 's/#.*//' -e 's/:.*//p' <origins.txt) \
	| sort >"$@"

${BUILD_ROOT}/all-files.txt: ${BUILD_ROOT}/all-projects.txt
	(builtin cd "${PROJECTS_DIR}"; \
	 find `cat ${BUILD_ROOT}/all-projects.txt` ! -path ".git*" -a -type f -print \
	 | sort) >"$@"

${BUILD_ROOT}/all-files.txt.gz: ${BUILD_ROOT}/all-files.txt

# Some files that match the fortran-file-patterns.txt list are not Fortran files.
# Eliminate them explicitly.
${BUILD_ROOT}/all-fortran-files.txt: fortran-file-patterns.txt \
	    fortran-exceptions.txt \
	    ${BUILD_ROOT}/all-files.txt
	grep -i -f fortran-file-patterns.txt ${BUILD_ROOT}/all-files.txt \
	| fgrep -v -f fortran-exceptions.txt >"$@"
	wc -l "$@"

# Compute the attributes in parallel and concatenate results.
# This is awkward, but cuts the time to generate them down dramatically.
# Note that since all-fortran-files.txt is sorted, the attributes
# file will be as well, without having to explicitly sort it.
# We touch it at the end, for the unusual circumstance
# where the file looks out of date, but no changes were
# made via `cpif,.
${BUILD_ROOT}/all-fortran-files-attr.txt:	${BUILD_ROOT}/all-fortran-files.txt \
	    ../bin/determine-attributes
	../bin/create-stats-file ${PROJECTS_DIR} \
	    ${BUILD_ROOT}/all-fortran-files.txt \
	    ../bin/determine-attributes \
	| ../bin/cpif "$@"
	wc -l "$@"

${BUILD_ROOT}/all-fortran-files-attr.txt.gz: ${BUILD_ROOT}/all-fortran-files-attr.txt

# All fixed-form Fortran files
${BUILD_ROOT}/all-fortran-files-fixed.txt:	${BUILD_ROOT}/all-fortran-files-attr.txt
	gawk  -F '\t' '/form:fixed/ { print $$1 }' \
	      <"$<" >"$@"
	wc -l "$@"

# All free-form Fortran files
${BUILD_ROOT}/all-fortran-files-free.txt:	${BUILD_ROOT}/all-fortran-files-attr.txt
	gawk -F '\t' '/form:free/ { print $$1 }' \
	      <"$<" >"$@"
	wc -l "$@"

# Results of running "file" on each Fortran file.
# This is useful for finding non-Fortran files that are
# disguised with Fortran-file-like names.
${BUILD_ROOT}/all-fortran-files-type.txt:	${BUILD_ROOT}/all-fortran-files.txt
	../bin/create-stats-file ${PROJECTS_DIR} \
	    ${BUILD_ROOT}/all-fortran-files.txt \
	    /usr/bin/file \
	| ../bin/cpif "$@"
	touch "$@"
	wc -l "$@"

${BUILD_ROOT}/all-fortran-files-type.txt.gz: ${BUILD_ROOT}/all-fortran-files-type.txt

# Line count of Fortran files in each project
${BUILD_ROOT}/all-fortran-files-lc.txt:	${BUILD_ROOT}/all-fortran-files-attr.txt
	gawk -F '\t' \
	    '{ file = $$1; \
	         lc = gensub(/.*lines:([0-9]*).*/, "\\1", 1, $$0)+0; \
	         printf("%8d\t%s\n", lc, file); \
	     }' \
	      <${BUILD_ROOT}/all-fortran-files-attr.txt >"$@"
	wc -l "$@"

${BUILD_ROOT}/all-projects-lc.txt:   ${BUILD_ROOT}/all-fortran-files-lc.txt
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
	      <${BUILD_ROOT}/all-fortran-files-lc.txt >"$@"
	wc -l "$@"

# List number of Fortran files in each project
${BUILD_ROOT}/all-projects-fortran-file-count.txt:	${BUILD_ROOT}/all-fortran-files-attr.txt
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
	      <${BUILD_ROOT}/all-fortran-files-attr.txt >"$@"
	wc -l "$@"

# Print some moderately interesting stats about the repositories.
${BUILD_DIR}/stats.txt:  ${BUILD_ROOT}/all-projects.txt \
	    ${BUILD_ROOT}/all-projects-lc.txt \
	    ${BUILD_ROOT}/all-files.txt \
	    ${BUILD_ROOT}/all-fortran-files.txt \
	    ${BUILD_ROOT}/all-fortran-files-lc.txt \
	    ${BUILD_ROOT}/all-projects-fortran-file-count.txt
	@(printf "%'12d projects\n" $$(wc -l <${BUILD_ROOT}/all-projects.txt); \
	  printf "%'12d files\n" $$(wc -l <${BUILD_ROOT}/all-files.txt); \
	  printf "%'12d Fortran files\n" $$(wc -l <${BUILD_ROOT}/all-fortran-files.txt); \
	  printf "%'12d Fortran lines\n" \
		$$(gawk -F '\t' \
	            '{ sum += gensub(/.*lines:([0-9]*).*/, "\\1", 1, $$0)+0 } \
                    END { print sum }' \
	            <${BUILD_ROOT}/all-fortran-files-attr.txt); \
	  printf "%'12d Fortran 77 lines\n" \
		$$(gawk -F '\t' \
	            '/form:fixed/ { \
                         sum += gensub(/.*lines:([0-9]*).*/, "\\1", 1, $$0)+0} \
                     END { print sum }' \
		<${BUILD_ROOT}/all-fortran-files-attr.txt); \
	  printf "%'12d cpp directive lines\n" \
	        $$(gawk -F '\t' \
	            '/cpreprocessor/ { \
                         sum += gensub(/.*cpreprocessor:([0-9]*).*/, "\\1", 1, $$0)+0 \
	             } \
                     END { print sum }' \
	           <${BUILD_ROOT}/all-fortran-files-attr.txt); \
	  echo "Largest projects:"; \
	  sort -rn -k 1 ${BUILD_ROOT}/all-projects-lc.txt \
		| head -10 \
		| awk $$'{ printf "%\'12d %s\\n", $$1, $$2}'; \
	  echo "Largest Fortran files:"; \
	  sort -rn -k 1 ${BUILD_ROOT}/all-fortran-files-lc.txt \
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

${BUILD_ROOT}/fortran-lang-projects.txt:	fortran-lang-category-urls.txt
	for P in `cat fortran-lang-category-urls.txt`; do \
		curl -L --no-progress-meter "$$P" \
		| sed -n -e '/^<p class="rubric".*<span class="/s;.*href="https*://\([^"]*\).*;\1;p'; \
	done \
	| sort -u >"$@"


# Get Fortran projects from https://fortran-lang.org/
# and add new ones as submodules of this project.
# These are organized into separate pages based on
# project category, so we gather the category pages first.

FORTRAN_LANG_URL = https://fortran-lang.org/packages
${BUILD_ROOT}/fortran-lang-category-urls.txt:	FORCE
	curl -L --no-progress-meter "${FORTRAN_LANG_URL}" 2>&1 \
	| sed -n -e '/^<h2><a class="reference internal" href="\([^"]*\).*/s;;${FORTRAN_LANG_URL}/\1;p' \
	| sort -u >"$@"

# Get new projects from github.com/Beliavsky/Fortran-code-on-GitHub.
# Only look at Git projects (not sourceforge or whatever).
beliavsky-new-projects:	${BUILD_ROOT}/beliavsky-projects.txt
	for P in `cat ${BUILD_ROOT}/beliavsky-projects.txt`; do \
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
${BUILD_ROOT}/beliavsky-projects.txt:	FORCE
	curl -L --no-progress-meter "${BELIAVSKY_URL}" \
		| sed -n -e 's;^<p dir="auto"><a href="https://\([^"]*\)".*;\1;p' \
	        | sed -e 's/#.*//' \
		| sort -u >"$@"

clean: FORCE
	rm -f ${BUILD_FILES}

FORCE:
