test:
	pytest --disable-pytest-warnings --cov=`pwd | rev | cut -d'/' -f1 | rev` ./tests
	rm -rf .coverage*

lint:
	pre-commit run --all-files

doc:
	tool=`pwd | rev | cut -d'/' -f1 | rev`
	rm -rf docs/source/$tool*
	rm -rf docs/source/auto_examples
	rm -rf docs/source/gen_modules
	cd docs && $(MAKE) clean html
	cd docs/build/html && python -m http.server 8080
