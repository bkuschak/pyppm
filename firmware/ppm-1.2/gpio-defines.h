
/* gpio-defines.h: PPM firmware general purpose input/output header.
 * Copyright (C) 2014  Bradley Worley  <geekysuavo@gmail.com>
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

/* ensure proper inclusion. */
#ifndef __PPM_GPIO_DEFINES_H__
#define __PPM_GPIO_DEFINES_H__

/* * * * GPIO HARDWARE CONNECTIONS * * * */

/* define connections for: light-emitting diodes. */
#define GPIO_LED_DDR  DDRC
#define GPIO_LED_PORT PORTC
#define GPIO_LED_PINA PORTC4
#define GPIO_LED_PINB PORTC5
#define GPIO_LED_PINC PORTC6

/* define connections for: constant current sink. */
#define GPIO_CCSEN_DDR  DDRC
#define GPIO_CCSEN_PORT PORTC
#define GPIO_CCSEN_PIN  PORTC2

/* define connections for: polarization relay. */
#define GPIO_RELAY_DDR  DDRD
#define GPIO_RELAY_PORT PORTD
#define GPIO_RELAY_PIN  PORTD3

/* define connections for: analog signal chain power. */
#define GPIO_LDOEN_DDR  DDRB
#define GPIO_LDOEN_PORT PORTB
#define GPIO_LDOEN_PIN  PORTB6

/* define connections for: hardware spi analog-to-digital converter. */
#define GPIO_HWSPI_DDR  DDRB
#define GPIO_HWSPI_SS   DDB0
#define GPIO_HWSPI_SCK  DDB1
#define GPIO_HWSPI_MOSI DDB2
/* GPIO_HWSPI_MISO need not be defined: master mode forces it as input. */

/* define connections for: software spi digital-to-analog converter. */
#define GPIO_SWSPI_DDR  DDRD
#define GPIO_SWSPI_PORT PORTD
#define GPIO_SWSPI_SS   PORTD0
#define GPIO_SWSPI_SCK  PORTD1
#define GPIO_SWSPI_MOSI PORTD2

/* define connections for: analog-to-digital converter chip-select. */
#define GPIO_ADCSEL_DDR  DDRD
#define GPIO_ADCSEL_PORT PORTD
#define GPIO_ADCSEL_PIN  PORTD6

/* * * * GPIO SOFTWARE FUNCTIONS: LEDs * * * */

/* gpio_led_init: initialize the state of the led outputs. */
#define gpio_led_init() \
  GPIO_LED_DDR |= ((1 << GPIO_LED_PINA) | \
                   (1 << GPIO_LED_PINB) | \
                   (1 << GPIO_LED_PINC)); \
  gpio_led_comm_off (); \
  gpio_led_adc_off (); \
  gpio_led_ccs_off ();

/* gpio_led_comm_on: turn on the first led output. */
#define gpio_led_comm_on() \
  GPIO_LED_PORT |= (1 << GPIO_LED_PINA);

/* gpio_led_comm_off: turn off the first led output. */
#define gpio_led_comm_off() \
  GPIO_LED_PORT &= ~(1 << GPIO_LED_PINA);

/* gpio_led_adc_on: turn on the second led output. */
#define gpio_led_adc_on() \
  GPIO_LED_PORT |= (1 << GPIO_LED_PINB);

/* gpio_led_adc_off: turn off the second led output. */
#define gpio_led_adc_off() \
  GPIO_LED_PORT &= ~(1 << GPIO_LED_PINB);

/* gpio_led_ccs_on: turn on the third led output. */
#define gpio_led_ccs_on() \
  GPIO_LED_PORT |= (1 << GPIO_LED_PINC);

/* gpio_led_ccs_off: turn off the third led output. */
#define gpio_led_ccs_off() \
  GPIO_LED_PORT &= ~(1 << GPIO_LED_PINC);

/* * * * GPIO SOFTWARE FUNCTIONS: CCS * * * */

/* gpio_ccs_init: initializes the state of the current source. */
#define gpio_ccs_init() \
  gpio_ccs_disable ();

/* gpio_ccs_enable: enable the polarization current source. this function
 * sends the CCS_EN pin into a high-Z state so it has no effect on the
 * voltage set-point at the CCS opamp non-inverting input.
 */
#define gpio_ccs_enable() \
  GPIO_CCSEN_DDR &= ~(1 << GPIO_CCSEN_PIN); \
  GPIO_CCSEN_PORT &= ~(1 << GPIO_CCSEN_PIN);

/* gpio_ccs_disable: disable the polarization current source. this function
 * ties the CCS_EN pin low to short the bottom leg of the CCS set-point
 * voltage divider to ground.
 */
#define gpio_ccs_disable() \
  GPIO_CCSEN_DDR |= (1 << GPIO_CCSEN_PIN); \
  GPIO_CCSEN_PORT &= ~(1 << GPIO_CCSEN_PIN);

/* * * * GPIO SOFTWARE FUNCTIONS: RELAY * * * */

/* gpio_relay_init: initializes the state of the POL/ACQ relay. */
#define gpio_relay_init() \
  GPIO_RELAY_DDR |= (1 << GPIO_RELAY_PIN); \
  gpio_relay_acq ();

/* gpio_relay_pol: engages the POL/ACQ relay into the POL mode by sending
 * the RELAY_EN pin high. this software requires that the default relay
 * position is ACQ mode.
 */
#define gpio_relay_pol() \
  GPIO_RELAY_PORT |= (1 << GPIO_RELAY_PIN);

/* gpio_relay_acq: disengages the POL/ACQ relay into the ACQ mode by sending
 * the RELAY_EN pin low. this software requires that the default relay
 * position is ACQ mode.
 */
#define gpio_relay_acq() \
  GPIO_RELAY_PORT &= ~(1 << GPIO_RELAY_PIN);

/* * * * GPIO SOFTWARE FUNCTIONS: LDO * * * */

/* gpio_ldo_init: initializes the state of the ASC power supply. */
#define gpio_ldo_init() \
  GPIO_LDOEN_DDR |= (1 << GPIO_LDOEN_PIN); \
  gpio_ldo_disable ();

/* gpio_ldo_enable: enables the analog signal chain power supply by
 * sending the LDO_EN pin high.
 */
#define gpio_ldo_enable() \
  GPIO_LDOEN_PORT |= (1 << GPIO_LDOEN_PIN);

/* gpio_ldo_disable: disables the analog signal chain power supply by
 * sending the LDO_EN pin low.
 */
#define gpio_ldo_disable() \
  GPIO_LDOEN_PORT &= ~(1 << GPIO_LDOEN_PIN);

/* * * * GPIO SOFTWARE FUNCTIONS: HW SPI * * * */

/* gpio_hwspi_init: initializes the hardware spi gpio ports. */
#define gpio_hwspi_init() \
  GPIO_HWSPI_DDR |= ((1 << GPIO_HWSPI_SS) | \
                     (1 << GPIO_HWSPI_SCK) | \
                     (1 << GPIO_HWSPI_MOSI));

/* * * * GPIO SOFTWARE FUNCTIONS: SW SPI * * * */

/* gpio_swspi_init: initializes the software spi gpio ports. */
#define gpio_swspi_init() \
  GPIO_SWSPI_DDR |= ((1 << GPIO_SWSPI_SS) | \
                     (1 << GPIO_SWSPI_SCK) | \
                     (1 << GPIO_SWSPI_MOSI)); \
  GPIO_SWSPI_PORT |= (1 << GPIO_SWSPI_SCK); \
  gpio_swspi_mosi_low ();

/* gpio_swspi_clock: toggles the clock pin of the software spi bus. */
#define gpio_swspi_clock() \
  GPIO_SWSPI_PORT &= ~(1 << GPIO_SWSPI_SCK); \
  GPIO_SWSPI_PORT |= (1 << GPIO_SWSPI_SCK);

/* gpio_swspi_mosi_low: sets the master output pin low. */
#define gpio_swspi_mosi_low() \
  GPIO_SWSPI_PORT &= ~(1 << GPIO_SWSPI_MOSI);

/* gpio_swspi_mosi_high: sets the master output pin high. */
#define gpio_swspi_mosi_high() \
  GPIO_SWSPI_PORT |= (1 << GPIO_SWSPI_MOSI);

/* * * * GPIO SOFTWARE FUNCTIONS: ADC /CS * * * */

/* gpio_adc_init: initialize the ADC /CS line state. */
#define gpio_adc_init() \
  gpio_hwspi_init (); \
  GPIO_ADCSEL_DDR |= (1 << GPIO_ADCSEL_PIN); \
  gpio_adc_deselect ();

/* gpio_adc_select: sends the ADC /CS line low to begin a conversion. */
#define gpio_adc_select() \
  GPIO_ADCSEL_PORT &= ~(1 << GPIO_ADCSEL_PIN);

/* gpio_adc_deselect: sends the ADC /CS line low to end a conversion. */
#define gpio_adc_deselect() \
  GPIO_ADCSEL_PORT |= (1 << GPIO_ADCSEL_PIN);

/* * * * GPIO SOFTWARE FUNCTIONS: DAC /CS * * * */

/* gpio_dac_init: initialize the DAC /CS line state. */
#define gpio_dac_init() \
  gpio_swspi_init (); \
  gpio_dac_deselect ();

/* gpio_dac_select: sends the DAC /CS line low to begin a transmission. */
#define gpio_dac_select() \
  GPIO_SWSPI_PORT &= ~(1 << GPIO_SWSPI_SS);

/* gpio_dac_deselect: sends the DAC /CS line high to end a transmission. */
#define gpio_dac_deselect() \
  GPIO_SWSPI_PORT |= (1 << GPIO_SWSPI_SS);

#endif /* __PPM_GPIO_DEFINES_H__ */
