PYTHON = venv/bin/python
GILLESPIE = $(PYTHON) -m src.gillespie.infrastructure.cli

# ── Moran (existing) ──────────────────────────────────────────────────────────
test:
	pytest tests

test-cov:
	pytest --cov=src --cov-report=html tests

run-deterministic:
	$(PYTHON) main.py --config configs/deterministic.yaml

run-deterministic-mutation:
	$(PYTHON) main.py --config configs/deterministic_mutation.yaml

run-bernoulli-mutation:
	$(PYTHON) main.py --config configs/bernoulli_mutation.yaml

run-all:
	$(MAKE) run-deterministic
	$(MAKE) run-deterministic-mutation
	$(MAKE) run-bernoulli-mutation

# ── Gillespie examples ────────────────────────────────────────────────────────
# Scenario 1: homeostasis — birth rate == death rate, no mutation.
# Population fluctuates around N0 and may go extinct.
gillespie-homeostasis:
	$(GILLESPIE) \
	  --N0 50 \
	  --lambda0 0.25 --mu0 0.25 --nu0 0.00 \
	  --T-max 20 --seed 42 \
	  --save-history results_gillespie_homeostasis.csv \
	  --save-debug results_gillespie_homeostasis_debug.csv

# Scenario 2: tumour growth with genomic instability and clonal evolution.
# lambda > mu drives net growth; nu0 > 0 generates mutant subclones that
# accumulate instability over time.
gillespie-tumour-growth:
	$(GILLESPIE) \
	  --N0 20 \
	  --lambda0 0.35 --mu0 0.20 --nu0 0.01 \
	  --instability-0 0.1 --buildup-0 0.0 \
	  --base-instability-buildup 0.005 \
	  --mutation-instability-jump 0.05 \
	  --fitness-gain 0.02 \
	  --T-max 30 --seed 7 \
	  --save-history results_gillespie_tumour_growth.csv \
	  --save-clones  results_gillespie_tumour_growth_clones.csv \
	  --save-debug results_gillespie_tumour_growth_debug.csv \
	  --top 15

# Scenario 3: logistic crowding — tumour growth with carrying capacity K0.
# Population grows until resource competition (crowding) limits expansion,
# modelling competition for space and nutrients inside a tumour niche.
gillespie-crowding:
	$(GILLESPIE) \
	  --N0 20 \
	  --lambda0 0.35 --mu0 0.20 --nu0 0.01 \
	  --use-logistic \
	  --K0 500 \
	  --instability-0 0.1 \
	  --base-instability-buildup 0.005 \
	  --mutation-instability-jump 0.05 \
	  --fitness-gain 0.02 \
	  --T-max 30 --seed 7 \
	  --save-history results_gillespie_crowding.csv \
	  --save-clones  results_gillespie_crowding_clones.csv \
	  --save-debug results_gillespie_crowding_debug.csv \
	  --top 15

gillespie-all:
	$(MAKE) gillespie-homeostasis
	$(MAKE) gillespie-tumour-growth
	$(MAKE) gillespie-crowding

# ── Environment ───────────────────────────────────────────────────────────────
create_venv:
	python -m venv venv

activate_venv:
	source venv/bin/activate

install-notebook:
	pip install -e .