test:
	pytest --cov=trisicell ./tests
	rm -rf .coverage*

lint:
	pre-commit run --all-files

doc:
	rm -rf docs/source/trisicell*
	rm -rf docs/source/auto_examples
	rm -rf docs/source/gen_modules
	cd docs && $(MAKE) clean html
	cd docs/build/html && python -m http.server 8080
