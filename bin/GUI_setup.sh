#setup
mamba create -n snakeui python=3.11 -y
mamba activate snakeui
mamba install -c conda-forge nodejs
pip install fastapi uvicorn websockets pydantic

curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/master/install.sh | bash

cd bin/frontend/

nvm install --lts
nvm use --lts

npm install express cors js-yaml
