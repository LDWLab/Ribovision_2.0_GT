# DESIRE

## Installation

1. Initial set up:

	a. Linux distribution (Native Linux system or <a href="https://docs.microsoft.com/en-us/windows/wsl/install-win10#update-to-wsl-2">WSL</a>)
	
	b. <a href=" https://code.visualstudio.com/docs/setup/setup-overview">Visual Studio Code</a> (VS Code)

	c. Follow the first two sections of <a href="https://code.visualstudio.com/docs/setup/setup-overview">this guide</a> to install a virtual python environment for the Django project (at least python3.6).

	d. Instal git and get an account on github.com. Contact the current admin to receive access to the project repository,
	as well as user account and password for the MySQL database.

	e. In the Linux startup file (.bashrc or .bash_profile) add these lines:

	>export DJANGO_SECRET_KEY='' (provided by admin)
	
	>export DJANGO_USERNAME='' (MySQL username provided by admin)
	
	>export DJANGO_PASSWORD='' (MySQL password provided by admin)

	f. Set up GaTech VPN with the Cisco AnyConnect Secure Mobility Client. Follow <a href="https://faq.oit.gatech.edu/content/how-do-i-get-started-campus-vpn">these</a> instructions.

2. Clone the <a href="https://github.com/LDWLab/DESIRE.git">project repository</a> in a new folder.

3. Using the command line from the root directory of the project rn the following commands:

	a. Activate the virtual environment

	>source ./env/bin/activate

	b. Install python requirements

	>pip install -r requirements.txt

	c. Install nodejs and its package manager (npm)

	>sudo apt-get install nodejs npm

	d. Install nodejs requirements (listed in package.json)

	>npm install

	e. Initial set up of nodejs scripts

	>./node_modules/.bin/webpack --config webpack.config.js

	f. (Optional) If developing a nodejs script

	>npm run watch

4. Open a VS Code folder in the root directory of the project.

5. The debugger should now work while connected to the GaTech VPN.