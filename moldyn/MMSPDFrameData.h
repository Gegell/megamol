/*
 * MMSPDFrameData.h
 *
 * Copyright (C) 2011 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOLCORE_MMSPDFRAMEDATA_H_INCLUDED
#define MEGAMOLCORE_MMSPDFRAMEDATA_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#include "vislib/Array.h"
#include "vislib/RawStorage.h"
#include "vislib/types.h"


namespace megamol {
namespace core {
namespace moldyn {


    /**
     * MMSPD file frame data
     */
    class MMSPDFrameData {
    public:

        /**
         * Class holding all particles of one type
         */
        class Particles {
        public:

            /**
             * Ctor
             */
            Particles(void);

            /**
             * Dtor
             */
            ~Particles(void);

            /**
             * Gets the number of particles
             *
             * @return The number of particles
             */
            inline UINT64 Count(void) {
                return this->count;
            }

            /**
             * Accesses the particle data
             *
             * @return The particle data
             */
            inline vislib::RawStorage& Data(void) {
                return this->data;
            }

            /**
             * Gets the particle data
             *
             * @return The particle data
             */
            inline const vislib::RawStorage& GetData(void) const {
                return this->data;
            }

            /**
             * Sets the number of particles
             *
             * @param cnt The number of particles
             */
            inline void SetCount(UINT64 cnt) {
                this->count = cnt;
            }

            /**
             * Test for equality
             *
             * @param rhs The right-hand side operand
             *
             * @return True if 'this' and 'rhs' are equal
             */
            bool operator==(const Particles& rhs);

        private:

            /** The number of particles */
            UINT64 count;

            /**
             * The particle data
             *
             * This is the raw data following the variable fields defined by
             * the particle type. (never storing the type id field!)
             */
            vislib::RawStorage data;

        };

        /**
         * Ctor
         *
         * @param owner The owning data module
         */
        MMSPDFrameData(void);

        /**
         * Dtor
         */
        virtual ~MMSPDFrameData(void);

        /**
         * Accesses the particle data
         *
         * @return The particle data
         */
        inline vislib::Array<Particles>& Data(void) {
            return this->data;
        }

        /**
         * Gets the particle data
         *
         * @return The particle data
         */
        inline const vislib::Array<Particles>& GetData(void) const {
            return this->data;
        }

        /**
         * Accesses the particle index reconstruction data
         *
         * @return The particle index reconstruction data
         */
        inline vislib::RawStorage& IndexReconstructionData(void) {
            return this->idxRec;
        }

        /**
         * Gets the particle index reconstruction data
         *
         * @return The particle index reconstruction data
         */
        inline const vislib::RawStorage& GetIndexReconstructionData(void) const {
            return this->idxRec;
        }

    private:

        /**
         * The particle data
         *
         * For each particle type (array index) the list of particles is stored
         * with their raw data corresponding to the field of the particle type
         * definition (optionally including UINT64 index; but never including
         * the type).
         */
        vislib::Array<Particles> data;

        /**
         * The particle index reconstruction data
         *
         * If loaded, this memory contains the blocks from the different arrays
         * of 'data' in the order they were stored in the data file. The first
         * integer defines the type index and the second integer defines the
         * number of particles.
         * All integers are RLE-Bit-Encoded.
         * (worst case 2 byte per particle)
         */
        vislib::RawStorage idxRec;

    };


} /* end namespace moldyn */
} /* end namespace core */
} /* end namespace megamol */

#endif /* MEGAMOLCORE_MMSPDFRAMEDATA_H_INCLUDED */
