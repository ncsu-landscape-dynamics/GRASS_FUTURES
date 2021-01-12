/*!
   \file keyvalue.c

   \brief Function to handle key-value pairs

   Derived from the Key_Value (lib/gis/key_value1.c) functions in GRASS
   GIS library which handled char*-char* pairs.

   (C) 2016 by Vaclav Petras and the GRASS Development Team

   This program is free software under the GNU General Public License
   (>=v2).  Read the file COPYING that comes with GRASS for details.

   \author Vaclav Petras
 */

#include <stdlib.h>
#include <string.h>
#include <grass/gis.h>

#include "keyvalue.h"

/*!
   \brief Allocate and initialize KeyValueIntInt structure

   \return poiter to allocated KeyValueIntInt structure
 */
struct KeyValueIntInt *KeyValueIntInt_create()
{
    struct KeyValueIntInt *kv = (struct KeyValueIntInt *) G_malloc(sizeof(struct KeyValueIntInt));
    G_zero(kv, sizeof(struct KeyValueIntInt));

    return kv;
}

/*!
   \brief Set value for given key

   \param[in,out] kv KeyValueIntInt structure to be modified
   \param key key to be set up
   \param value value for given key
 */
void KeyValueIntInt_set(struct KeyValueIntInt *kv, int key, int value)
{
    int n;

    for (n = 0; n < kv->nitems; n++)
        if (key == kv->key[n])
            break;

    if (n == kv->nitems) {
        if (n >= kv->nalloc) {
            size_t size;

            if (kv->nalloc <= 0)
                kv->nalloc = 8;
            else
                kv->nalloc *= 2;

            size = kv->nalloc * sizeof(int);
            kv->key = (int *) G_realloc(kv->key, size);
            kv->value = (int *) G_realloc(kv->value, size);
        }

        kv->key[n] = key;
        kv->value[n] = value;
        kv->nitems++;
        return;
    }

    kv->value[n] = value;
}

/*!
   \brief Find given key

   \param key key to be found
   \param[out] value pointer to store the value if found
   \param kv pointer to KeyValueIntInt structure

   \returns TRUE if the key is found (and sets the value)
   \returns FALSE if no key found
 */
int KeyValueIntInt_find(const struct KeyValueIntInt *kv, int key, int *value)
{
    int n;

    if (!kv)
        return FALSE;

    for (n = 0; n < kv->nitems; n++)
        if (key == kv->key[n]) {
            *value = kv->value[n];
            return TRUE;
        }

    return FALSE;
}

/*!
   \brief Free allocated KeyValueIntInt structure

   \param[in] kv KeyValueIntInt structure to be G_freed
 */
void KeyValueIntInt_free(struct KeyValueIntInt *kv)
{
    if (!kv)
        return;

    G_free(kv->key);
    G_free(kv->value);
    kv->nitems = 0;                /* just for safe measure */
    kv->nalloc = 0;
    G_free(kv);
}

/*!
   \brief Allocate and initialize KeyValueIntFloat structure

   \return poiter to allocated KeyValueIntFloat structure
 */
struct KeyValueIntFloat *KeyValueIntFloat_create()
{
    struct KeyValueIntFloat *kv = (struct KeyValueIntFloat *) G_malloc(sizeof(struct KeyValueIntFloat));
    G_zero(kv, sizeof(struct KeyValueIntFloat));

    return kv;
}

/*!
   \brief Set value for given key

   \param[in,out] kv KeyValueIntFloat structure to be modified
   \param key key to be set up
   \param value value for given key
 */
void KeyValueIntFloat_set(struct KeyValueIntFloat *kv, int key, float value)
{
    int n;

    for (n = 0; n < kv->nitems; n++)
        if (key == kv->key[n])
            break;

    if (n == kv->nitems) {
        if (n >= kv->nalloc) {
            size_t size;

            if (kv->nalloc <= 0)
                kv->nalloc = 8;
            else
                kv->nalloc *= 2;

            size = kv->nalloc * sizeof(int);
            kv->key = (int *) G_realloc(kv->key, size);
            kv->value = (float *) G_realloc(kv->value, kv->nalloc * sizeof(float));
        }

        kv->key[n] = key;
        kv->value[n] = value;
        kv->nitems++;
        return;
    }

    kv->value[n] = value;
}

/*!
   \brief Find given key

   \param key key to be found
   \param[out] value pointer to store the value if found
   \param kv pointer to KeyValueIntFloat structure

   \returns TRUE if the key is found (and sets the value)
   \returns FALSE if no key found
 */
int KeyValueIntFloat_find(const struct KeyValueIntFloat *kv, int key, float *value)
{
    int n;

    if (!kv)
        return FALSE;

    for (n = 0; n < kv->nitems; n++)
        if (key == kv->key[n]) {
            *value = kv->value[n];
            return TRUE;
        }

    return FALSE;
}

/*!
   \brief Free allocated KeyValueIntFloat structure

   \param[in] kv KeyValueIntFloat structure to be G_freed
 */
void KeyValueIntFloat_free(struct KeyValueIntFloat *kv)
{
    if (!kv)
        return;

    G_free(kv->key);
    G_free(kv->value);
    kv->nitems = 0;                /* just for safe measure */
    kv->nalloc = 0;
    G_free(kv);
}