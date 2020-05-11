#ifndef FUTURES_KEYVALUE_H
#define FUTURES_KEYVALUE_H

struct KeyValueIntInt
{
    int nitems;
    int nalloc;
    int *key;
    int *value;
};

struct KeyValueCharInt
{
    int nitems;
    int nalloc;
    const char **key;
    int *value;
};

struct KeyValueIntInt *KeyValueIntInt_create();
void KeyValueIntInt_set(struct KeyValueIntInt *kv, int key, int value);
int KeyValueIntInt_find(const struct KeyValueIntInt *kv, int key, int *value);
void KeyValueIntInt_free(struct KeyValueIntInt *kv);


struct KeyValueCharInt *KeyValueCharInt_create();
void KeyValueCharInt_set(struct KeyValueCharInt *kv, const char *key, int value);
int KeyValueCharInt_find(const struct KeyValueCharInt *kv, const char *key, int *value);
void KeyValueCharInt_free(struct KeyValueCharInt *kv);
#endif /* FUTURES_KEYVALUE_H */
