notebook:
	pip3 install jupyter
	jupyter notebook --port 8001 --allow-root --ip=0.0.0.0 ""

docker:
	docker build --cache-from rvdata:latest --tag rvdata:latest .
	docker run -it rvdata:latest bash

regression_tests:
	pytest -x --cov=core --cov=instruments --pyargs tests.regression
	coveralls

.PHONY: init
