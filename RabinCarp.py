from Bio import SeqIO

# функция для вычисления хеша строки
def hash_string(s):
    p = 31 # простое число
    m = 10**9 + 9 # большое простое число
    h = 0
    for c in s:
        h = (h * p + ord(c)) % m
    return h

# функция для вычисления хешей всех подстрок длины L в строке s
def hash_substrings(s, L):
    p = 31 # простое число
    m = 10**9 + 9 # большое простое число
    h = 0
    pL = pow(p, L, m)
    for i in range(L):
        h = (h * p + ord(s[i])) % m
    hashes = [h]
    for i in range(L, len(s)):
        h = (h * p - ord(s[i-L]) * pL + ord(s[i])) % m
        hashes.append(h)
    return hashes

# функция для поиска подстроки P в строке T с помощью алгоритма Рабина-Карпа
def rabin_karp(p, t):
    n = len(t) - len(p) + 1
    p_hash = hash_string(p)
    t_hashes = hash_substrings(t, len(p))
    for i in range(n):
        if p_hash == t_hashes[i]:
            if p == t[i:i+len(p)]:
                return i+1
    return -1

# чтение данных из файлов
p_record = SeqIO.read('P.fa', 'fasta')
t_record = SeqIO.read('T.fa', 'fasta')

# поиск вхождения P в T
locus = rabin_karp(str(p_record.seq), str(t_record.seq))
print('Alignment locus:', locus)