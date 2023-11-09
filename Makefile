all:	all-projects.txt all-files.txt all-fortran-files.txt \
			all-projects-lc.txt all-fortran-files-lc.txt \
			all-projects-fortran-file-count.txt stats.txt

update-all: update all

# Update projects from their remote repos
update:
	git submodule update --remote

# Create a list of all Fortran projects we have,
# all files, and all Fortran files.
all-projects.txt:	Makefile origins.txt .gitmodules
	(git submodule foreach -q 'echo $$sm_path'; \
	 sed -n -e 's/#.*//' -e 's/:.*//p' <origins.txt) \
    | sort >"$@"

all-files.txt:	Makefile all-projects.txt
	find `cat all-projects.txt` ! -path ".git*" -a -type f -print \
	| sort >"$@"

FIND_IS_FORTRAN= -type f \( -iname "*.f" -o -iname "*.f[0-9][0-9]" \
	-o -iname "*.for" -o -name "*.fpp"" \
	-o -iname "*.ftn" -o -iname "*.fypp \)

all-fortran-files.txt:	all-projects.txt Makefile
	find `cat all-projects.txt` ${FIND_IS_FORTRAN} -print \
	| sort >"$@"

# Word count each Fortran file
all-fortran-files-lc.txt:	all-fortran-files.txt Makefile
	tr '\n' '\0' <all-fortran-files.txt \
	| xargs -0 wc -l \
	| grep -v 'total$$' >"$@"

# Word count each Fortran file
all-projects-lc.txt:	all-fortran-files.txt Makefile
	for P in `cat all-projects.txt`; do \
		echo "$$(find $$P ${FIND_IS_FORTRAN} -print0 \
			| xargs -0 cat \
			| wc -l)\t$$P"; \
	done >"$@"

# List number of Fortran files in each project
all-projects-fortran-file-count.txt:	all-fortran-files.txt Makefile
	for P in `cat all-projects.txt`; do \
		echo "$$P\t$$(find $$P ${FIND_IS_FORTRAN} -print \
			| wc -l | tr -d ' ')"; \
	done >"$@"

# Print some moderately interesting stats about the repositories.
stats.txt:  all-projects.txt all-projects-lc.txt all-files.txt \
			all-fortran-files.txt all-fortran-files-lc.txt \
			all-projects-fortran-file-count.txt
	@(printf "%'12d projects\n" $$(wc -l <all-projects.txt); \
	  printf "%'12d files\n" $$(wc -l <all-files.txt); \
	  printf "%'12d Fortran files\n" $$(wc -l <all-fortran-files.txt); \
	  printf "%'12d Fortran lines\n" \
			$$(awk '{ sum += $$1} END { print sum}' <all-fortran-files-lc.txt); \
	  printf "%'12d Fortran 77 lines\n" \
			$$(awk '/\.[Ff]$$/ { sum += $$1} END { print sum}' \
			<all-fortran-files-lc.txt); \
	  printf "%'12d cpp directive lines\n" \
	        $$(tr '\n' '\0' <all-fortran-files.txt \
			   | xargs -0 cat \
	           | egrep '^[[:space:]]*#' | wc -l); \
	  echo "Largest projects:"; \
	  sort -rn -k 1 all-projects-lc.txt \
		| head -10 \
		| awk $$'{ printf "%\'12d %s\\n", $$1, $$2}'; \
	  echo "Largest Fortran files:"; \
	  sort -rn -k 1 all-fortran-files-lc.txt \
		| head -10 \
		| awk $$'{ printf "%\'12d %s\\n", $$1, $$2}'; \
	  echo "Projects with no Fortran files:"; \
	  grep '[^0-9]0$$' all-projects-fortran-file-count.txt \
	    | sed -e 's/^/  /') \
	 | tee "$@"


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
	    if grep "^$$D:" exceptions.txt >/dev/null; then \
			echo "$$D is on the exception list."; \
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
	    if grep "^$$D:" exceptions.txt >/dev/null; then \
			echo "$$D is on the exception list."; \
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

BELIAVSKY_URL = https://github.com/Beliavsky/Fortran-code-on-GitHub
beliavsky-projects.txt:	FORCE
	curl -L --no-progress-meter "${BELIAVSKY_URL}" \
		| sed -n -e 's;^<p dir="auto"><a href="https://\([^"]*\)".*;\1;p' \
		| sort -u >"$@"


FORCE:
