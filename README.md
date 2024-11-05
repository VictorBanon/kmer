## How to install

At the root of this repository, run:
```
# Set up virtual environment and activate it
python3 -m venv venv
source venv/bin/activate
```

```
pip install maturin
cd rust
maturin develop --release
cd ..
# Install package
pip install .
run
```