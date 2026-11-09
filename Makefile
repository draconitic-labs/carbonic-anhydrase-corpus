.PHONY: build refresh offline validate test

build:
	.venv/bin/python scripts/build_corpus.py

refresh:
	.venv/bin/python scripts/build_corpus.py --refresh

offline:
	.venv/bin/python scripts/build_corpus.py --offline

validate:
	.venv/bin/python scripts/validate_corpus.py

test: validate
	.venv/bin/python -m unittest discover -s tests -v

