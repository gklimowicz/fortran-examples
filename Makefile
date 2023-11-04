FORTRAN_LANG_URL = https://fortran-lang.org/en/packages
FIND_IS_FORTRAN= -type f \( -iname "*.f" -o -iname "*.f[0-9][0-9]" -o -iname "*.ftn" \)

all:	all-projects.txt all-files.txt all-fortran-files.txt \
			all-projects-lc.txt all-fortran-files-lc.txt \
			stats.txt

update-all: update all

# Update projects from their remote repos
update:
	git submodule update --remote

# Create a list of all Fortran projects we have,
# all files, and all Fortran files.
all-projects.txt:	Makefile origins.txt .gitmodules
	(git submodule foreach -q 'echo $$sm_path'; \
	 sed -n -e 's/#.*//' -e 's/:.*//p' <origins.txt) \
    | sort > "$@"

all-files.txt:	Makefile all-projects.txt
	find `cat all-projects.txt` ! -path ".git*" -a -type f -print \
	| sort >"$@"

all-fortran-files.txt:	all-projects.txt Makefile
	find `cat all-projects.txt` ${FIND_IS_FORTRAN} -print \
	| sort >"$@"

# Word count each Fortran file
all-fortran-files-lc.txt:	all-fortran-files.txt Makefile
	cat all-fortran-files.txt \
	| tr '\n' '\0' \
	| xargs -0 wc -l \
	| grep -v 'total$$' >"$@"

# Word count each Fortran file
all-projects-lc.txt:	all-fortran-files.txt Makefile
	for P in `cat all-projects.txt`; do \
		echo "$$(find $$P ${FIND_IS_FORTRAN} -print0 \
			| xargs -0 cat \
			| wc -l)\t$$P"; \
	done >"$@"

# Print some moderately interesting stats about the repositories.
stats.txt:  all-projects.txt all-projects-lc.txt all-files.txt \
			all-fortran-files.txt all-fortran-files-lc.txt
	@(printf "%'12d projects\n" $$(wc -l <all-projects.txt); \
	  printf "%'12d files\n" $$(wc -l <all-files.txt); \
	  printf "%'12d Fortran files\n" $$(wc -l <all-fortran-files.txt); \
	  printf "%'12d Fortran lines\n" \
			$$(awk '{ sum += $$1} END { print sum}' <all-fortran-files-lc.txt); \
	  printf "%'12d Fortran 77 lines\n" \
			$$(awk '/\.[Ff]$$/ { sum += $$1} END { print sum}' \
			<all-fortran-files-lc.txt); \
	  printf "%'12d cpp directive lines\n" \
	        $$(tr '\n' '\0' <all-fortran-files.txt | xargs -0 cat | grep '^#' | wc -l); \
	  echo "Largest projects:"; \
	  sort -rn -k 1 all-projects-lc.txt \
		| head -10 \
		| awk $$'{ printf "%\'12d %s\\n", $$1, $$2}'; \
	  echo "Largest Fortran files:"; \
	  sort -rn -k 1 all-fortran-files-lc.txt \
		| head -10 \
		| awk $$'{ printf "%\'12d %s\\n", $$1, $$2}') \
	 | tee "$@"


# Add new projects we may find lying about.
add-new-projects:	fortran-lang-new-projects

# Get new projects for fortran-lang.org.
# Only look at Git projects (not sourceforge or whatever).
fortran-lang-new-projects:	fortran-lang-projects
	for P in `cat fortran-lang-projects`; do \
		case "$$P" in \
			git*) ;; \
			*) continue;; \
		esac; \
		D=$${P/*\/}; \
		if [[ -d "$$D" ]]; then \
			: echo "$$D exists already"; \
		else \
			echo "$$D new"; \
			git submodule add ssh://git@$$P \
			|| git submodule add https://$$P; \
			git submodule update --init "$$D"; \
		fi; \
	done

fortran-lang-projects.txt:	fortran-lang-category-urls.txt
	for P in `cat fortran-lang-category-urls.txt`; do \
		curl -L --no-progress-meter "$$P" \
		| sed -n -e '/^<p class="rubric"><span class="/s;.*href="https*://\([^"]*\).*;\1;p'; \
	done \
	| sort -u >"$@"

fortran-lang-category-urls.txt:	FORCE
	curl -L --no-progress-meter "${FORTRAN_LANG_URL}" 2>&1 \
	| sed -n -e '/^<h2><a class="reference internal" href="\([^"]*\).*/s;;${FORTRAN_LANG_URL}/\1;p' \
	| sort -u >"$@"

FORCE:
