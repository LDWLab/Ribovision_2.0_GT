/**
* Copyright 2018, Plotly, Inc.
* All rights reserved.
*
* This source code is licensed under the MIT license found in the
* LICENSE file in the root directory of this source tree.
*/

import { createProvider } from 'react-redux';

import { storeKey } from './storeOptions';

export default createProvider(storeKey);
