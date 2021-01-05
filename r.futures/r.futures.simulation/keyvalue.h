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

struct KeyValueIntFloat
{
    int nitems;
    int nalloc;
    int *key;
    float *value;
};

struct KeyValueIntInt *KeyValueIntInt_create();
void KeyValueIntInt_set(struct KeyValueIntInt *kv, int key, int value);
int KeyValueIntInt_find(const struct KeyValueIntInt *kv, int key, int *value);
void KeyValueIntInt_free(struct KeyValueIntInt *kv);


struct KeyValueCharInt *KeyValueCharInt_create();
void KeyValueCharInt_set(struct KeyValueCharInt *kv, const char *key, int value);
int KeyValueCharInt_find(const struct KeyValueCharInt *kv, const char *key, int *value);
void KeyValueCharInt_free(struct KeyValueCharInt *kv);

struct KeyValueIntFloat *KeyValueIntFloat_create();
void KeyValueIntFloat_set(struct KeyValueIntFloat *kv, int key, float value);
int KeyValueIntFloat_find(const struct KeyValueIntFloat *kv, int key, float *value);
void KeyValueIntFloat_free(struct KeyValueIntFloat *kv);
#endif /* FUTURES_KEYVALUE_H */
