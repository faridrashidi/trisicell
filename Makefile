test:
	pytest --disable-pytest-warnings --cov=`pwd | rev | cut -d'/' -f1 | rev` ./tests
	rm -rf .coverage*

lint:
	pre-commit run --all-files

docs:
	tool=`pwd | rev | cut -d'/' -f1 | rev`
	cd docs
	rm -rf ./source/$tool*
	make clean html
	cd build/html
	serve
	cd ../../..
