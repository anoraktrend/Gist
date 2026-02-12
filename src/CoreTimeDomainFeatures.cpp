//=======================================================================
/** @file CoreTimeDomainFeatures.cpp
 *  @brief Implementations of common time domain audio features
 *  @author Adam Stark
 *  @copyright Copyright (C) 2013  Adam Stark
 *
 * This file is part of the 'Gist' audio analysis library
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
//=======================================================================

#include "CoreTimeDomainFeatures.h"
#include <cmath>

//===========================================================
template <class T>
CoreTimeDomainFeatures<T>::CoreTimeDomainFeatures()
{
}

//===========================================================
template <class T>
T CoreTimeDomainFeatures<T>::rootMeanSquare (const std::vector<T>& buffer)
{
    T sum = 0;
    const size_t n = buffer.size();
    // Use x * x instead of pow(x, 2) for better performance
    for (size_t i = 0; i < n; i++)
    {
        T sample = buffer[i];
        sum += sample * sample;
    }
    return std::sqrt (sum / static_cast<T>(n));
}

//===========================================================
template <class T>
T CoreTimeDomainFeatures<T>::peakEnergy (const std::vector<T>& buffer)
{
    T peak = -10000.0;
    const size_t n = buffer.size();
    for (size_t i = 0; i < n; i++)
    {
        T absSample = std::fabs(buffer[i]);
        if (absSample > peak)
            peak = absSample;
    }
    return peak;
}

//===========================================================
template <class T>
T CoreTimeDomainFeatures<T>::zeroCrossingRate (const std::vector<T>& buffer)
{
    T zcr = 0;
    const size_t n = buffer.size();
    // Optimized: check sign change more efficiently
    for (size_t i = 1; i < n; i++)
    {
        // Check if signs differ (one positive, one negative or zero)
        if ((buffer[i] > 0) != (buffer[i - 1] > 0))
        {
            zcr += 1.0;
        }
    }
    return zcr;
}

//===========================================================
template class CoreTimeDomainFeatures<float>;
template class CoreTimeDomainFeatures<double>;
