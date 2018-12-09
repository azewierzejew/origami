(* Autor: Antoni Żewierżejew *)
(* Reviewer: Robert Michna *)


(** Punkt na płaszczyźnie *)
type point = float * float


(** Poskładana kartka: ile razy kartkę przebije szpilka wbita w danym punkcie *)
type kartka = point -> int


(** mały epsilon, żeby uwzględniać niedokładność floata *)
let epsilon = 1e-9


(** [nalezy v1 v2 v] zwraca czy [v] należy do przedziały (v1, v2) (co do epsilona) *)
let nalezy v1 v2 v = (v1 -. epsilon <= v) && (v <= v2 +. epsilon) 


(** [prostokat p1 p2] zwraca kartkę, reprezentującą domknięty prostokąt
o bokach równoległych do osi układu współrzędnych i lewym dolnym rogu [p1]
a prawym górnym [p2]. Punkt [p1] musi więc być nieostro na lewo i w dół
od punktu [p2]. Gdy w kartkę tę wbije się szpilkę wewnątrz
(lub na krawędziach) prostokąta, kartka zostanie przebita 1 raz,
w pozostałych przypadkach 0 razy *)
let prostokat (x1, y1) (x2, y2) = fun (x, y) -> if nalezy x1 x2 x && nalezy y1 y2 y then 1 else 0


(** [odlegosc p1 p2] zwraca odleglość pomiędzy tymi dwoma punktami *)
let odleglosc (x1, y1) (x2, y2) = hypot (x2 -. x1) (y2 -. y1)


(** [pole p1 p2 p3] zwraca dwukrotność pola trójkąta o takich wierzchołkach *) 
let pole (x1, y1) (x2, y2) (x3, y3) = 
    x1 *. y2 +. x2 *. y3 +. x3 *. y1 -. (x1 *. y3 +. x2 *. y1 +. x3 *. y2)


(** [kolko p r] zwraca kartkę, będącą kółkiem domkniętym o środku w punkcie [p] i promieniu [r] *)
let kolko p r = fun pnt ->
    if r < 0. then raise (Invalid_argument "Negative radius in [kolko]") else
    if odleglosc p pnt <= r +. epsilon then 1 else 0

(** [odbicie p1 p2 p] zwraca odbicie punktu [p] względem prostej zawierającej [p1] i [p2] *)
let odbicie ((x1, y1) as p1) ((x2, y2) as p2) (xp, yp) = 
    let a = y1 -. y2 in
    let b = x2 -. x1 in
    let c = x1 *. y2 -. x2 *. y1 in
    let odl = (a *. xp +. b *. yp +. c) /. (hypot a b) in
    let dl = odleglosc p1 p2 in
    let (xw, yw) = (a *. odl *. (-2.) /. dl, b *. odl *. (-2.) /. dl) in
    (xp +. xw, yp +. yw)


(** [zloz p1 p2 k] składa kartkę [k] wzdłuż prostej przechodzącej przez
punkty [p1] i [p2] (muszą to być różne punkty). Papier jest składany
w ten sposób, że z prawej strony prostej (patrząc w kierunku od [p1] do [p2])
jest przekładany na lewą. Wynikiem funkcji jest złożona kartka. Jej
przebicie po prawej stronie prostej powinno więc zwrócić 0.
Przebicie dokładnie na prostej powinno zwrócić tyle samo,
co przebicie kartki przed złożeniem. Po stronie lewej -
tyle co przed złożeniem plus przebicie rozłożonej kartki w punkcie,
który nałożył się na punkt przebicia. *)
let zloz p1 p2 kar = fun p ->
    if odleglosc p1 p2 < epsilon then raise (Invalid_argument "Points too close in [zloz]") else
    let wyz = pole p1 p2 p in
    if abs_float (wyz /. (odleglosc p1 p2)) < epsilon then kar p else
    if wyz < 0. then 0 else
    (kar p) + (kar (odbicie p1 p2 p))

(** [skladaj [(p1_1,p2_1);...;(p1_n,p2_n)] k = zloz p1_n p2_n (zloz ... (zloz p1_1 p2_1 k)...)]
czyli wynikiem jest złożenie kartki [k] kolejno wzdłuż wszystkich prostych z listy *) 
let skladaj l kar = List.fold_left (fun k (p1, p2) -> zloz p1 p2 k) kar l
