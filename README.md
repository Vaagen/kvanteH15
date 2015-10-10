# kvanteH15
Numeric solution to schrodinger equation




GIT COMMANDS

git clone https://github.com/Vaagen/kvanteH15.git
	Kloner mappen fra github

git add <filnavn> 
	Legger til endringer i <filnavn>

git add --all <mappenavn>
	legger til mappe med alle filene

git commit -m ‘Forklaring’ 
	gjør klar endringene og lagrer dem sammen med forklaringen

git pull origin
	Denne laster inn endringer i origin

git status
		viser hva som er lagt til og fjernet av filer>

NB: Bruk pull før du bruker push, da lagrer du ikke over noe av det som ligger på origin…

git push origin	
	Laster opp endringene til origin

git log --oneline
	Viser liste over tidligere endringer

git checkout <commit number>
	tar deg til tilstanden ved denne commit-en

git checkout master
	tar deg tilbake til nyeste tilstand

git fetch
git reset --hard origin
	Disse to i kombinasjon fjerner det du har og tar over Origin slik den er nå

git clean 
    Fjerner alle untracked files
    -d  for å ta med untracked files i undermapper
    -n  for å få oversikt over hva som vil slettes hvis du kjører git clean
    -f  for å tvinge git til å slette filene

git show-branch
    Viser info om branches, se https://wincent.com/wiki/Understanding_the_output_of_%22git_show-branch%22