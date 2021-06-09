test:
	pytest --disable-pytest-warnings --cov=`pwd | rev | cut -d'/' -f1 | rev` ./tests
	rm -rf .coverage*

lint:
	pre-commit run --all-files

doc:
	tool=`pwd | rev | cut -d'/' -f1 | rev`
	cd docs
	rm -rf ./source/$tool*
	rm -rf ./source/auto_examples
	rm -rf ./source/gen_modules
	make clean html
	cd build/html
	serve
	cd ../../..
