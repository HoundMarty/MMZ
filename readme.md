# Otázky ke zkoužce


## Biologické a akustické signály - úvod a opakování

### Teoretické

1. [Jaký je původ biosignálů a co je akční potenciál buňky?](otazky\1T1.md)
2. [Uveďte příklady pravidelných a nepravidelných biosignálů.](otazky\1T2.md)
3. [Které technické parametry biosignálů a měřící technicky jsou pro zpracování určující?](otazky\1T3.md)
4. [Jaké třídy filtrů jsou realizovatelné?](otazky\1T4.md)
5. [Které filtry mají (zobecněnou) lineární fázi a za jakých okolností tuto vlastnost filtru požadujeme?](otazky\1T5.md)

### Praktické:

1. [Odstranění plovoucí stejnosměrné složky v signálu pomocí filtru. Popište princip. Jaký filtr je vhodné zvolit?](otazky\1P1.md)
2. [Odstranění harmonických složek pomocí filtru. Popište postup. Jaký filtr je vhodné zvolit?](otazky\1P2.md)


## Kvadratická kritéria pro porovnávání signálů

### Teoretické

1. Náhodná veličina: definice a vlastnosti a její role v modelování náhodných signálů.
2. Definice a odhad střední hodnoty a rozptylu signálu.
3. Korelace a kovariance signálů a jejich odhad. 
4. Korelační koeficient: definice a vlastnosti.
5. Jaká je souvislost mezi kvadratickou vzdáleností signálů, jejich rozptyly a korelací?
6. Křížová korelace: definice a vlastnosti.

### Praktické

1. Metoda synchronizovaného průměrování: účel, princip a realizace.


## Optimální filtry ve smyslu kvadratické vzdálenosti: Wienerův filtr

### Teoretické

1. Least Mean Square (LMS) nebo též Wienerův filtr: princip a odvození.
2. Jaká je role autokovariance a křížové korelace ve výpočtu LMS filtru?

### Praktické

1. Odhad vzájemného zpoždění signálů pomocí LMS filtru: princip a realizace.
2. Potlačení (akustického) echa pomocí LMS filtru: popis problému, princip řešení a realizace.


## Optimální filtry ve smyslu kvadratické vzdálenosti: Adaptivní filtry

### Teoretické

1. Adaptivní LMS filtr: definice problému, princip a odvození.
2. Adaptivní RLS algoritmus: princip.
3. Použití adaptivních filtrů ve frekvenční oblasti: princip.

### Praktické

1. Odhad vzájemného zpoždění signálů pomocí adaptivních filtrů.


## Komprimované vzorkování

### Teoretické

1. Co je řídký signál a jak definujeme kritérium řídkosti? Co je oblast řídkosti signálu? Je každý signál v nějaké oblasti řídký?
2. Popište princip ztrátové komprese a rekonstrukce signálu v řídké oblasti. 
3. Použití adaptivních filtrů ve frekvenční oblasti: princip.
4. Komprimované vzorkování: princip.
5. Koherence bází: definice. Je pro komprimované vzorkování signálu lepší jsou-li měřící a řídká báze málo nebo hodně koherentní? Proč?
6. Co říká teorém o minimálním počtu vzorkům, které je třeba změřit, aby byla možná úplná rekonstrukce signálu? Jaké jsou funkční závislosti na délce signálu, řídkosti a koherenci bází?

### Praktické

1. Rekonstrukce signálu pomocí Orthogonal Matching Pursuit (OMP) nebo l1-magic algoritmu.


## Signály z více senzorů a jejich zpracování

### Teoretické

1. Lineární modely šíření signálů v prostoru: definice a účel a popis ve frekvenční oblasti.
2. Vícekanálový filtr: definice a popis v časové a frekvenční oblasti.
3. Přenosová funkce vícekanálového filtru v prostředí bez odrazů.
4. Delay-and-sum beamformer: definice a vlastnosti.
5. MVDR beamformer: definice a princip.
6. Optimální vícekanálový filtr ve smyslu kvadratické vzdálenosti: princip a definice.

### Praktické

1. Výpočet vícekanálového LMS filtru v časové oblasti.


## Analýza hlavních komponent

### Teoretické

1. Komponenta vícekanálového signálu: definice, rozptyl komponenty a souvislost s kovarianční maticí vícekanálového signálu.
2. Směrový vektor signálu: definice.
3. Hlavní komponenty vícekanálového signálu: definice a výpočet.
4. Vlastnosti hlavních komponent vícekanálového signálu.
5. Redukce dimenze vícekanálového signálu pomocí PCA: princip a postup.

### Praktické

1. Použití PCA: dekompozice a rekonstrukce vícekanálového signálu.


## Slepá separace signálů

### Teoretické

1. Definice úlohy slepé separace signálů; okamžitý model míšení signálů; přeurčený, určený a nedourčený model.
2. Princip analýzy nezávislých komponent, nejednoznačnost řešení.
3. Vzájemná informace: definice a vlastnosti.

### Praktické

1. Použití ICA: dekompozice a rekonstrukce vícekanálového signálu. 
