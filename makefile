test:
	pytest tests

test-cov:
	pytest --cov=src --cov-report=html tests

run-deterministic:
	python main.py --config configs/deterministic.yaml

run-deterministic-mutation:
	python main.py --config configs/deterministic_mutation.yaml

run-bernoulli-mutation:
	python main.py --config configs/bernoulli_mutation.yaml

run-all:
	$(MAKE) run-deterministic
	$(MAKE) run-deterministic-mutation
	$(MAKE) run-bernoulli-mutation

activate_venv:
	source venv/bin/activate